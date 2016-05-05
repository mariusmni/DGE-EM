//TO KEEP THE COMPILER HAPPY:
// import scala.collection.JavaConversions._
// import edu.uconn.engr.dna.util.ImplicitScalaJavaConversions._

package edu.uconn.engr.dna.dgesim;
import edu.uconn.engr.dna.util.IsoformCleavageSites

import edu.uconn.engr.dna.sim.Read
import edu.uconn.engr.dna.util.Utils._
import edu.uconn.engr.dna.random._
import edu.uconn.engr.dna.format._
import edu.uconn.engr.dna.sim.io.DefaultOutputWriter
import edu.uconn.engr.dna.sim.format.{ FixedNamingConvention, FASTQReadFormatter }
import java.io.FileInputStream
import edu.uconn.engr.dna.io.{ DefaultTwoFieldParser, GTFParser }
import edu.uconn.engr.dna.probability.ProbabilityDistribution

import edu.uconn.engr.dna.util.ImplicitScalaJavaConversions._
import scala.collection.JavaConversions._
import edu.uconn.engr.dna.util.{ DefaultRandomAccessMap, StringToDoubleRandomAccessMap, RandomAccessMap , EmUtils}

import scala.collection.breakOut

/**
 * User: marius
 * Date: Jul 26, 2010
 * Time: 2:14:49 PM
 *
 */
object DGESimStartup {
	def println(m: String) {
		Console.err.println("sim:" + m)
	}

	def main(args: Array[String]) {
		var simulationPropertiesFile = "simulation.properties"
		if (args.length == 1) {
			simulationPropertiesFile = args(0)
		} else if (args.length > 1) {
			println("Arguments: [simulationPropertiesFile]")
			println("If no argument is specified, the file \"simulation.properties\" will be used")
			return
		}
		run(simulationPropertiesFile)
	}

	def run(simulationPropertiesFile: String) {
		val start = System.currentTimeMillis

		println("Reading simulation properties...")
		val params = new SimulationParams(new FileInputStream(simulationPropertiesFile))
		println(params.toString)
		val prefix = createFileNamePrefix(params)

		println("Reading GTF file... ")
		val p: Tuple2[Clusters, Isoforms] = new GTFParser().parse(new FileInputStream(params.getGTF))
		var (clusters, isoforms) = p
		if (params.getClusters != null) {
			println("Parsing clusters file...");
			clusters = new Clusters;
			DefaultTwoFieldParser
			.getInvertedTwoFieldParser(clusters)
			.parse(new FileInputStream(params.getClusters));
			EmUtils.synchronizeClustersWithIsoforms(clusters, isoforms)
		};
		println("read " + isoforms.size + " isoforms and " + clusters.size + " clusters")

		println("Loading cluster abundances(weights)...")
		val clusterAbundances = loadClusterAbundances(params.clusterDistributionFileName)
		for (c <- clusters.idIterator)
			if (clusterAbundances.get(c) == null) {
				clusterAbundances.put(c, 0.0);
				println(c)
			}

		val tagLength = params.tagLength
		println("Reading known gene sequences...")
		// Isoforms sequences are loaded as they are in the file,
		// which is in coding strand orientation!
		val isoformSequences = new IsoformSequencesFromFile(params.isoformSequencesFile,
															isoforms, false)

		// add polyA
		val polyA = Array.fill(250)('A').toString
		for (tag <- isoformSequences.getAllTags.toList) {
			isoformSequences.put(tag, isoformSequences.getSequence(tag) + polyA);
		}

		IsoformCleavageSites.checkIsoforms(isoformSequences, isoforms, clusters)
		val sitesForIsoform: Map[String, IndexedSeq[Int]] = IsoformCleavageSites.
		computeIsoformCleaveSites(isoformSequences, params.cleavePatterns, tagLength)
		val aeSitePicker = sitesForIsoform.map(pair => 
			(pair._1, new GeometricProbabilityIndexPicker(params.p, params.randomNumberGeneratorSeed, pair._2.length)))
    
		val clusterPicker = new MapRandomPicker(new CustomWeightIdPicker[java.lang.String, Clusters](clusterAbundances,
																									 params.randomNumberGeneratorSeed),
												clusters)

		val isoformDistribution = params.isoformDistribution
		val isoformPicker = createIsoformPicker(isoformDistribution,
												isoforms, clusters, params.randomNumberGeneratorSeed)


//		val nc = if (params.printToConsole) null else new FixedNamingConvention(prefix + ".fastr")
//		val ow = new DefaultOutputWriter(
//			new FastrReadFormatter(isoformSequences),
//			nc, null, false)

		val nc = if (params.printToConsole) null else new FixedNamingConvention(prefix + ".fastq")
		val ow = new DefaultOutputWriter(
			new FASTQReadFormatter(isoformSequences, false),
			nc, null, false)

		var tagCount = params.numberOfTags
		println("generating " + tagCount + " tags...")
		var counter = 0
		var iterations = 0
		var N1 = 0
		while (counter < tagCount) {
			val cluster = clusterPicker.pick(clusters)
			if (cluster.size == 0) {
				println("WARNING: Cluster " + cluster.getName + " is empty")
			} else {
				val isoform = isoformPicker.pick(cluster)
				if (isoform != null) {
					val sites = sitesForIsoform(isoform.getName)
					if (sites.size > 0) {
						val siteIndex = aeSitePicker(isoform.getName).pick(null)
						if (siteIndex >= 0) {
							val site = sites(siteIndex)
							// polyA tailed sequences: tag always fits
//							if (tagLength + site <= isoform.length) {
							if (siteIndex == 0) {
								N1+=1
							}
							val tagName = "tag" + counter + "_" + isoform.getName + "_" + siteIndex
							ow.writeReads(new Read(tagName, cluster.getName, isoform, 1 + site, tagLength, isoform.getStrand))
							counter += 1
//							}
						}
					}
				}
			}
			iterations+=1
		}
		println("ACTUAL PICKS " + iterations)
		println("Picks from rightmost site" + N1)
		ow.close

		if (params.outputIsoformWeights) {
			println("writing isoform weights used in the simulation...")
			var globalIsoWeights = createGlobalIsoformWeights(isoformDistribution, clusters, clusterAbundances, isoforms)
			globalIsoWeights = sort(normalize(globalIsoWeights), false)
			writeValues(globalIsoWeights.entrySet, prefix + ".iso_weights")
			println("writing gene weights used in the simulation...")
			var globalGeneWeights = createGlobalIsoformWeights(isoformDistribution, clusters, clusterAbundances, isoforms)
			globalIsoWeights = sort(normalize(globalIsoWeights), false)
			writeValues(globalIsoWeights.entrySet, prefix + ".iso_weights")
		}
		println(String.format("done. (%.2fs)", ((System.currentTimeMillis - start) / 1000.0).asInstanceOf[AnyRef]))
	}

	private def createFileNamePrefix(params: SimulationParams) = "nt=" + params.numberOfTags +
	"_tl=" + params.tagLength + "_isofD=" + params.isoformDistribution + "_p=" + params.p

	private def createIsoformPicker(
		isoformDistribution: ProbabilityDistribution,
		isoforms: Isoforms,
		clusters: Clusters,
		randomNumberGeneratorSeed: Int): RandomPicker[Isoform, Cluster] = {
		val isoformWeights = new DefaultRandomAccessMap[String]
		for (cluster <- clusters.groupIterator) {
			var i: Int = 0
			for (isoform <- select(isoforms, cluster)) {
				isoformWeights.add(isoform.getName, isoformDistribution.getWeight(i, cluster.size))
				i += 1
			}
		}
		val weightedIDPicker = new CustomWeightIdPicker[String, Cluster](isoformWeights, randomNumberGeneratorSeed)
		new MapRandomPicker(weightedIDPicker, isoforms)
	}

	private def loadClusterAbundances(fileName: String): RandomAccessMap[java.lang.String, java.lang.Double] = {
		DefaultTwoFieldParser.getRegularTwoFieldParser(new StringToDoubleRandomAccessMap[String])
		.parse(new FileInputStream(fileName))
	}

	def createGlobalIsoformWeights(isoformDistribution: ProbabilityDistribution,
								   clusters: Clusters,
								   clusterAbundances: RandomAccessMap[java.lang.String, java.lang.Double],
								   isoforms: Isoforms): RandomAccessMap[java.lang.String, java.lang.Double] = {
		val isoWeights = new DefaultRandomAccessMap[java.lang.String]
		for (cluster <- clusters.groupIterator) {
			val clusterAbundance: Double = clusterAbundances.getValue(cluster.getName).doubleValue
			var i = 0
			for (isoform <- select(isoforms, cluster)) {
				val isoWeight = clusterAbundance * isoformDistribution.getWeight(i, cluster.size)
				isoWeights.add(isoform.getName, isoWeight)
				i += 1
			}
		}
		isoWeights
	}

}