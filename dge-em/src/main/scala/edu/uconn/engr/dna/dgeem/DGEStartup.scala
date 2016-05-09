package edu.uconn.engr.dna.dgeem

import java.util.concurrent.ArrayBlockingQueue
import edu.uconn.engr.dna.util.BatchThreadPoolExecutor
import edu.uconn.engr.dna.dgeem.DGEOptionParser._
import edu.uconn.engr.dna.format.{IsoformSequencesFromFile, MemoryMapChromosomeSequences, TaggedSequences, Clusters,
								  Isoforms}
import edu.uconn.engr.dna.io.{DefaultTwoFieldParser, GTFParser}
import edu.uconn.engr.dna.util._

import java.io.{InputStreamReader, FileReader, BufferedReader, FileInputStream}
import collection.mutable.ListBuffer
import edu.uconn.engr.dna.util.F
import collection.JavaConverters._



import joptsimple.OptionSet


import java.{lang => jl, util => ju}

/**
 * User: marius
 * Date: Jul 30, 2010
 */
object DGEStartup {
	private val BUFF_SIZE = 1 << 16;
	val DEBUG = false

	def debug(msg: ()=>String) {
		if (DEBUG) {
			println(msg())
		}
	}
	
	def main(args: Array[String]) {
		val options = parseAndValidateOptions(args);
		val (clusters, isoforms) = parseGTF(options);

		var genome: TaggedSequences = null;
		val timer = new TimingTool;
		timer.start(null);
		
		println("reading known gene sequences...");
		// Isoforms sequences are loaded as they are in the file,
		// which is in coding strand orientation!
		val isoformSequencesFile = options.valueOf(OP_MRNA).toString;
		val isoformSequences = new IsoformSequencesFromFile(isoformSequencesFile, isoforms, false);
		IsoformCleavageSites.checkIsoforms(isoformSequences, isoforms, clusters);
		val isoformToClusterMap = Utils.createIsoformToClusterMap(isoforms, clusters);

		val cleavePatterns = options.valueOf(OP_ENZIME).toString.split(",").toList;
		val runUniq = options.has(OP_UNIQ);

		var inputFiles =  options.nonOptionArguments.asInstanceOf[java.util.List[String]].asScala.toList;
		if (inputFiles.isEmpty) {
			inputFiles = List[String](null)
		}
		val nThreads = Runtime.getRuntime.availableProcessors;
		if (options.has(OP_PREPEND_ENZYME)) {
			if (cleavePatterns.size != 1) {
				println("ERROR: option --" + OP_PREPEND_ENZYME + " can only be used for single enzymes");
				System.exit(1);
			}
		}
		val mism = if (options.has(OP_MISMATCHES))
			options.valueOf(OP_MISMATCHES).asInstanceOf[jl.Integer].intValue
		else
			1;
		val tagMatcher = new TagMatcher(isoformSequences, cleavePatterns,
										mismatches=mism, qualScores=true,
										enzymeIncluded=(!options.has(OP_PREPEND_ENZYME)))
//			val tagMatcher = new IsoformHash(isoformSequences, cleavePatterns)
//			val factory = F.factory(() => new DGEParserParameterRunnable(tagMatcher/*, classes*/))
//				val cluster = new UnionFindClustering[WeightedTagClassJ](
//					F.converter(tagClass => jc.asJavaIterator(tagClass.isoforms.iterator)))
		val cls = new ReadClassCompactAndClusterParameterRunnable[WeightedTagClassJ](
			F.converter(tagClass => tagClass.isoforms.toIterator.asJava.asInstanceOf[ju.Iterator[jl.Object]]),
			F.binaryOperator((toBeRemoved, toBeKept) => {
					toBeKept.multiplicity += toBeRemoved.multiplicity;
					null.asInstanceOf[jl.Void]}))
		val cluster = if (runUniq) {
			new ParameterRunnable[jl.Iterable[WeightedTagClassJ], ju.List[ju.List[WeightedTagClassJ]]] {
				override def run(x: jl.Iterable[WeightedTagClassJ]) {
					cls.run(keepUniq(x, isoformToClusterMap));
				}
				override def done = cls.done
			}
		} else {
			cls
		}
		val clusterProcess = new SingleBatchThreadPoolExecutor(cluster, math.max(100, 10 * nThreads));
		if (options.has(OP_LIMIT_NTAGS)) {
			val limit = options.valueOf(OP_LIMIT_NTAGS).asInstanceOf[jl.Integer].intValue
			println("!Limited to max " + limit  + " tags")
			DGEParseAndForward.limitNTags =  limit
		}
//		val tagCollector = new TagCollector
		val factory = F.factory(() => new DGEParseAndForward(tagMatcher,
															 clusterProcess
//															 tagCollector
			))
		for (file <- inputFiles) {
			try {
				if (DGEParseAndForward.notFull)  {
					if (file == null)
						println("Parsing tags from standard input...")
					else
						println("Parsing tags from " + file + "... ")

					val inputReader = new BufferedReader(
						if (file != null) new java.io.FileReader(file)
						else new InputStreamReader(System .in), BUFF_SIZE);
													 
					if (file == null || getExtension(file).toLowerCase == ".fastq" || getExtension(file).toLowerCase == ".fq") {
					  new ThreadPooledFastqParser(nThreads, 0x1000, 0x10000, factory).parse(inputReader)
					} else {
						new ThreadPooledFastrParser(nThreads, 0x1000, 0x10000, factory).parse(inputReader)
					}
				}
			} catch {
				case e: Exception => {
						e.printStackTrace(System.out)
					}
			}
		}
		
		if (tagMatcher.noMatch > 0) {
			println("! Could not find a match for " + tagMatcher.noMatch + " tags")
		}
		val nrCleaveSites = tagMatcher.nrCleaveSites
		val isoformsHit = nrCleaveSites.count(p => p._2 > 0)
		debug(()=>"Isoforms with at least one site " + isoformsHit)
		debug(()=>"Total sites for patterns " + (nrCleaveSites.foldRight(0) {
					(p, count) => count + p._2
				}));

		val nrCleaveSitesJavaMap = new ju.HashMap[jl.String, jl.Integer]();
		for ((k, v) <- nrCleaveSites) {
			nrCleaveSitesJavaMap.put(k, v);
		}

		val tagClusters = clusterProcess.waitForTermination;
		debug(()=>"Done clustering. Clusters " + tagClusters.size);
		debug(()=>"Largest cluster " + tagClusters.asScala.foldLeft(0)((m, l) => math.max(m, l.size)))
		debug(()=>"Tags in clusters " + tagClusters.asScala
				.foldLeft(0)((c, l) => c + l.size))
		println(()=>"Tags: " + tagClusters.asScala
				.foldLeft(0)((c, l) => c + l.asScala.foldLeft(0)((s, t) => s + t.multiplicity)))

		var tags = Utils.flatten(tagClusters)
		if (runUniq) {
			debug(()=>"Uniq tag classes " + tags.size)
			println("Uniq tags: " + tags.asScala.foldLeft(0)((s, t) => s + t.multiplicity))
		}
		var freq: ju.Map[jl.String, jl.Double] = if (runUniq) {
			val map = new ju.HashMap[jl.String, jl.Double]
			for (t <- tags.asScala;
				 gene = isoformToClusterMap.get(t.isoforms(0))) {
				val count = if (map.containsKey(gene)) map.get(gene).doubleValue
				else 0.0
				map.put(gene, count + t.multiplicity)
			}
			map
		} else {
			val em = new DGEEMParameterRunnable(isoforms, nrCleaveSitesJavaMap)
			em.run(tags)
			em.done
		}
		//printf("%d components did not converge\n", DGEEMJ.noConverge)

		var fileName = if (options.has(OP_OUTPUT_PREFIX))
			options.valueOf(OP_OUTPUT_PREFIX).toString + ".fastq"
		else if (inputFiles.size == 1) inputFiles.head
		else "dge.fastq";

		if (!runUniq) {
			// Uniq does not output isoform expression
			val isoOutputFileName = getOutputFilePrefix(fileName) + ".iso_estimates";
			println("Writing isoform FPKMs to " + isoOutputFileName);
			EmUtils.addMissingIds(freq, isoforms.idIterator);
			Utils.writeValues(Utils.sortEntriesDesc(freq), isoOutputFileName);
		}

		val geneOutputFileName = getOutputFilePrefix(fileName) + ".gene_estimates";
		println("Writing gene FPKMs to " + geneOutputFileName);
		if (!runUniq) {
			// Uniq already reported gene frequencies
			freq = Utils.groupByCluster(freq, isoformToClusterMap);
		}
		EmUtils.addMissingIds(freq, clusters.idIterator);
//		freq = EmUtils.normalize(freq);
		Utils.writeValues(Utils.sortEntriesDesc(freq), geneOutputFileName);
		/*
		 val collectedTags = tagCollector.output
		 println("Collected " + collectedTags.size + " tags")
		 val thresholds = List(1, 2, 4, 8, 16, 32, 64, 128, 256).map(_*1000000);
		 var previousThreshold = 0
		 for (threshold <- thresholds) {
		 if (previousThreshold < collectedTags.size) {
		 cluster.run(collectedTags.subList(previousThreshold, math.min(threshold, collectedTags.size)))
		 for (i <- previousThreshold until math.min(threshold, collectedTags.size)) {
		 collectedTags.set(i, null); // let gc collect
		 }
		 var collapsedTags = cluster.collectAllReadClasses();
		 //				println("Clustering " + classes.size + " classes ")
		 //				println("Total tags " + classes.foldLeft(0)((c, tc) => c + tc.multiplicity))

		 //				val tagClusters = clusterProcess.waitForTermination;
		 //				println("Done clustering. Clusters " + tagClusters.size);
		 //				println("Largest cluster " + tagClusters.asScala.foldLeft(0)((m, l) => math.max(m, l.size)))
		 //				println("Tags in clusters " + tagClusters.asScala
		 //						.foldLeft(0)((c, l) => c + l.size))
		 //				println("Tag multiplicity in clusters " + tagClusters.asScala
		 //						.foldLeft(0)((c, l) => c + l.asScala.foldLeft(0)((s, t) => s + t.multiplicity)))


		 println("Collapsed tags " + collapsedTags.size)
		 println("Collapsed tag multiplicity " + collapsedTags.asScala.
		 foldLeft(0)((s, t) => s + t.multiplicity))
		 println("Running EM")

		 var em = new DGEEMParameterRunnable(p, isoforms, nrCleaveSitesJavaMap)
		 //				var tags = Utils.flatten(tagClusters)
		 //				if (options.has(OP_UNIQ)) {
		 //					tags = keepUniq(tags, isoformToClusterMap)
		 //					println("Uniq tags " + tags.size)
		 //					println("Uniq tag multiplicity " + tags.asScala.foldLeft(0)((s, t) => s + t.multiplicity))
		 //				}
		 //				em.run(tags)
		 //				em.run(uniq)
		 em.run(collapsedTags)
		 var freq: ju.Map[jl.String, jl.Double] = em.done
		 //				printf("%d components did not converge\n", DGEEMJ.noConverge)

		 //				val maps = BatchThreadPoolExecutor
		 //				.newInstance(nThreads,
		 //							 new ArrayBlockingQueue[ju.List[WeightedTagClassJ]](math.max(100, 10 * nThreads)),
		 //							 F.factory(()=>new DGEEMParameterRunnable(p, isoforms, nrCleaveSitesJavaMap)),
		 //							 false)
		 //				.processAll(tagClusters)
		 //				.waitForTermination();


		 //				println("Collecting all frequencies")
		 //				var freq: ju.Map[jl.String, jl.Double] = new ju.HashMap[jl.String, jl.Double]();
		 //				val it = maps.iterator;
		 //				while (it.hasNext) {
		 //					val map = it.next;
		 //					for ((k, v) <- map) {
		 //						freq.put(k, v);
		 //					}
		 //				}


		 var fileName = if (options.has(OP_OUTPUT_PREFIX))
		 options.valueOf(OP_OUTPUT_PREFIX).toString + ".fastq"
		 else if (inputFiles.size == 1) inputFiles.head
		 else "dge.fastq";
		 fileName = fileName.replace(".fastq", "_" + threshold + ".fastq")
		 writeFreqs(fileName, freq, isoformToClusterMap, isoforms, clusters)

		 collapsedTags = keepUniq(collapsedTags, isoformToClusterMap)
		 println("Uniq tags " + collapsedTags.size)
		 println("Uniq tag multiplicity " + collapsedTags.asScala.
		 foldLeft(0)((s, t) => s + t.multiplicity))
		 println("Running EM Uniq")
		 em = new DGEEMParameterRunnable(p, isoforms, nrCleaveSitesJavaMap)
		 em.run(collapsedTags)
		 freq = em.done
		 writeFreqs(fileName.replace(".fastq", "uniq.fastq"), freq,
		 isoformToClusterMap, isoforms, clusters)
		 }
		 previousThreshold = threshold
		 }
		 */
		timer.stop;
		printf("Done. (%.2fs)\n", timer.getGlobalTime / 1000.0);
	};

	def keepUniq(tags: jl.Iterable[WeightedTagClassJ], iso2Cluster: ju.Map[jl.String, jl.String]) = {
		val filtered = new ju.ArrayList[WeightedTagClassJ]();
		for (t <- tags.asScala) {
			if (isUniq(t, iso2Cluster)) {
				filtered.add(t);
			}
		}
		filtered
	}

	def isUniq(t: WeightedTagClassJ, iso2Cluster: ju.Map[jl.String, jl.String])= {
		var uniq = true
		var cluster: jl.String = null
		for (i <- 0 until t.isoforms.length) {
			val c = iso2Cluster.get(t.isoforms(i))
			if (cluster == null) {
				cluster = c
			} else {
				if (!c.equals(cluster)) {
					uniq = false
				}
			}
		}
		uniq
	}

	def parseGTF(options: OptionSet): Tuple2[Clusters, Isoforms] = {
		println("Parsing GTF file...");
		val giParser = new GTFParser;
		val p = giParser.parse(new FileInputStream(options.valueOf(OP_GTF).toString));
		var clusters = p.getFirst;
		val isoforms = p.getSecond;
		if (options.has(OP_CLUSTERS_FILE)) {
			println("Parsing clusters file...");
			clusters = new Clusters;
			DefaultTwoFieldParser
			.getInvertedTwoFieldParser(clusters)
			.parse(new FileInputStream(options.valueOf(OP_CLUSTERS_FILE).toString));
			EmUtils.synchronizeClustersWithIsoforms(clusters, isoforms)
		};
		println("Found " + isoforms.size + " isoforms and " + clusters.size + " genes");
		(clusters, isoforms)
	}

	def getExtension(file: String) = file.substring(file.length - 6)

	def dropExtension(file: String) = file.replace(getExtension(file), "")

	def getOutputFilePrefix(file: String) = if (file != null) dropExtension(file) else "dge"

	def parseAndValidateOptions(args: Array[String]) = {
		val parser = new DGEOptionParser;
		val options: OptionSet = try {
			parser.parse(args: _*)
		} catch {
			case e: Exception => {
					e.printStackTrace;
					parser.printHelpOn(System.out);
					System.exit(1);
					null
				}
		};
		if (!options.has(OP_GTF)
			|| options.has(OP_HELP)
			|| !options.has(OP_MRNA)) {
			parser.printHelpOn(System.out);
			System.exit(1);
		};
		options
	}

};

object DGEParseAndForward {
	var limitNTags = Int.MaxValue
	var currentCount = 0
	val lock = new Object

	def notFull = {
		lock.synchronized {
			currentCount < limitNTags
		}
	}
}
class DGEParseAndForward(tagMatcher: TagMatcher,
						 forwardProcess: ParameterRunnable[jl.Iterable[WeightedTagClassJ], _])
extends ParameterRunnable[ju.List[String], Unit] {
	import DGEParseAndForward._
	def run(item: ju.List[String]) {

		// for "perfect" mapping
		//								val n = item.size
		//								for (i <- 0 until n by 2) {
		//									val tagName = item.get(i)
		//									val tag = item.get(i+1)
		//									val a = tagName.split('_')
		//									val isoform = a(1)
		//									val pos = a(2).toInt
		//									val tc = new TagClass(1, List((isoform, pos, 1.0)))
		//									tclasses.synchronized {
		//										tclasses.add(tc)
		//										if (tagLen == -1) {
		//											tagLen = tag.length
		//										}
		//									}
		//								}

//		println("Batch of " + item.size)
		if (notFull)  {
			val tcs = new ju.ArrayList[WeightedTagClassJ]
			val iterator = item.iterator
			while (iterator.hasNext) {
				val tag = iterator.next
				val qual = iterator.next
				val tc = tagMatcher.getTagClass(tag, qual)
				if (tc != null) {
					lock.synchronized {
						if (notFull) {
							tcs.add(tc)
							currentCount += 1
						}
					}
				}
			}
			if (!tcs.isEmpty) {
				forwardProcess.run(tcs)
			}
		}
	};

	def done() {}
}

class DGEEMParameterRunnable(isoforms: Isoforms, nrCleaveSites: ju.Map[jl.String, jl.Integer])
extends ParameterRunnable[ju.Collection[WeightedTagClassJ], ju.Map[jl.String, jl.Double]] {
	val freq = new ju.HashMap[jl.String, jl.Double]();

	def run(tags: ju.Collection[WeightedTagClassJ]) {
		val em = new DGEEMJ(isoforms, nrCleaveSites);
		val n = tags.size
		val array = new Array[WeightedTagClassJ](n)
		var i = 0
		var it = tags.iterator
		while (it.hasNext) {
			array(i) = it.next
			i += 1
		}
		freq.putAll(em.computeFreq(array));
	};

	def done() = freq; //jc.asScalaMap(freq).toMap;
}

