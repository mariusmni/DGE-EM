/*
 *  Copyright 2011 marius.
 * 
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 * 
 *       http://www.apache.org/licenses/LICENSE-2.0
 * 
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 *  under the License.
 */

package edu.uconn.engr.dna.dgeem

import edu.uconn.engr.dna.format.Clusters
import edu.uconn.engr.dna.format.IsoformSequencesFromFile
import edu.uconn.engr.dna.format.Isoforms
import edu.uconn.engr.dna.format.TaggedSequences
import edu.uconn.engr.dna.io.DefaultTwoFieldParser
import edu.uconn.engr.dna.io.GTFParser
import edu.uconn.engr.dna.util.EmUtils
import edu.uconn.engr.dna.util.IsoformCleavageSites
import java.io.FileInputStream
import java.io.FileReader
import scala.collection.JavaConverters._

object ComputeCleaveSites {

	def main(args: Array[String]) {
		val genes = new java.util.Properties;
		genes.load(new FileReader("qpcr_genes.txt"));
		//val (clusters, isoforms) = parseGTF("hg19KnownGene.gtf", "hg19KnownToGnfAtlas2.txt");
		val (clusters, isoforms) = parseGTF("hg19ENSGene.gtf", "hg19ENStranscriptToGene.txt");

		println("reading known gene sequences...");
		// Isoforms sequences are loaded as they are in the file,
		// which is in coding strand orientation!
		val isoformSequencesFile = "hg19ENSKnownGeneMrna-polyA200.txt";
		val isoformSequences = new IsoformSequencesFromFile(isoformSequencesFile, isoforms, false);
		IsoformCleavageSites.checkIsoforms(isoformSequences, isoforms, clusters);

//		val enzymes = Array("AATT", "ASST", "CGCG", "GGCC", "AGCT", "RGCY",
//							"TGCA", "YATR", "ACGT", "CATG", "CCGC,GCGG", "CCGG",
//							"CTAG", "GATC", "GCGC", "GTAC", "TCGA", "TTAA");
		val enzymes = Array("GATC", "CATG", "RGCY");
		for (i <- 0 until enzymes.length) {
			printPatterns(enzymes(i), isoformSequences, 21, clusters, //null)
//						  jc.asScalaIterable(genes.keySet).map(g => g.toString)
						List("ENSG00000181222"))
//			for (j <- (i+1) until enzymes.length) {
//				printPatterns(enzymes(i)+","+enzymes(j), isoformSequences, 21)
//			}
		}
	}

	def printPatterns(cleavePatterns: String,
					  isoformSequences: TaggedSequences,
					  tagLen: Int,
					  clusters: Clusters,
					  interestingGenes: Iterable[String]) {
		val nrCleaveSites = IsoformCleavageSites.sitesForIsoforms(
			isoformSequences, cleavePatterns.split(",").toList, tagLen);
		val isosCut = nrCleaveSites.count(p => p._2 > 0)

		if (interestingGenes != null) {
			println("Gene subset: " + interestingGenes.size)
			val interestingClusters = interestingGenes.map(g => clusters.get(g));
			for (g <- interestingGenes) {
				val cluster = clusters.get(g)
				for (i <- cluster.idIterator().asScala) {
					interestingGenes.foreach(g => println(i + " " + nrCleaveSites(i)))
				}
			}
			val genesMinOneIsoCut = interestingClusters.filter(c =>
				c.idIterator.asScala.exists(i => nrCleaveSites(i) > 0))

			val genesAllIsosCut = interestingClusters.filter(c =>
				c.idIterator.asScala.forall(i => nrCleaveSites(i) > 0))
			for (g <- genesAllIsosCut) {
				println(g)
			}
			println(cleavePatterns + "\t" + isosCut + "\t" +
				genesMinOneIsoCut.size + "\t" + genesAllIsosCut.size)
		} else {
			val genesMinOneIsoCut = clusters.groupIterator.asScala.count(c =>
				c.idIterator.asScala.exists(i => nrCleaveSites(i) > 0))
			val genesAllIsosCut = clusters.groupIterator.asScala.count(c =>
				c.idIterator.asScala.forall(i => nrCleaveSites(i) > 0))
			println(cleavePatterns + "\t" + isosCut + "\t" +
				genesMinOneIsoCut + "\t" + genesAllIsosCut)
		}
	}

	def parseGTF(gtf: String, clust: String): Tuple2[Clusters, Isoforms] = {
		println("Parsing GTF file...");
		val giParser = new GTFParser;
		val p = giParser.parse(new FileInputStream(gtf));
		var clusters = p.getFirst;
		val isoforms = p.getSecond;
		if (clust != null) {
			println("Parsing clusters file...");
			clusters = new Clusters();
			DefaultTwoFieldParser.getInvertedTwoFieldParser(clusters).parse(
							new FileInputStream(clust));
			EmUtils.synchronizeClustersWithIsoforms(clusters, isoforms);
		}
		println("Found " + isoforms.size + " isoforms and " + clusters.size + " genes");
		(clusters, isoforms)
	}
}
