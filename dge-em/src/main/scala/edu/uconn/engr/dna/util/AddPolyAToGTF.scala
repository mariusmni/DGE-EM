package edu.uconn.engr.dna.util

import edu.uconn.engr.dna.io.GTFParser
import java.io.FileInputStream
import java.io.PrintWriter
import scala.collection.JavaConversions._

object AddPolyAToGTF {
	def main(args: Array[String]) {
//		val args = Array("250", "knownGeneGnfAtlas2.gtf", "knownGeneGnfAtlas2-polyA250.gtf", 	"hg18_knownGene_GnfAtlas2")
		if (args.length != 4) {
			println("Arguments: tailLength gtfFile outputFile source")
			System.exit(-1)
		}
		val tailLen = args(0).toInt
		val gtfFile = args(1)
		val outputFile = args(2)
		val source = args(3)
		val r = new GTFParser().parse(new FileInputStream(gtfFile))
    val clusters = r.getFirst()
    val isoforms = r.getSecond()
		val isoformToClusterMap = Utils.createIsoformToClusterMap(isoforms, clusters)
		val maxCoord = isoforms.isoformIterator.map(i => i.getExons.getEnd).max
		val dummyExonStart = maxCoord + 10
		val dummyExonEnd = dummyExonStart + tailLen - 1
		val out = new PrintWriter(outputFile)
		val printFormat = "%s\t%s\texon\t%s\t%s\t0.000000\t%s\t.\tgene_id \"%s\"; transcript_id \"%s\";\n"
		for (isoform <- isoforms.isoformIterator;
				 cluster = isoformToClusterMap.get(isoform.getName);
				 exons = isoform.getExons) {

			if (isoform.getStrand == '-') {
				out.printf(printFormat, isoform.getChromosome, source,
									 (-tailLen).toString, "-1",
									 isoform.getStrand.toString, cluster, isoform.getName)
			}

			for (i <- 0 until exons.size;
					 start = exons.getStart(i);
					 end = exons.getEnd(i)) {
				out.printf(printFormat, isoform.getChromosome, source,
									 start.toString, end.toString, isoform.getStrand.toString,
									 cluster, isoform.getName)
			}

			if (isoform.getStrand == '+') {
				out.printf(printFormat, isoform.getChromosome, source,
									 dummyExonStart.toString, dummyExonEnd.toString,
									 isoform.getStrand.toString, cluster, isoform.getName)
			}
		}
		out.close
	}
}


