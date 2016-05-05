package edu.uconn.engr.dna.util

import edu.uconn.engr.dna.format.IsoformSequencesFromFile
import java.io.PrintWriter
import scala.collection.JavaConversions._
import scala.math._

object AddPolyAToMRNA {
	def main(args: Array[String]) {
		if (args.length != 4 && args.length != 3) {
			println("Arguments: tailLength mRnaFile outputFile [basesPerLine]")
			System.exit(-1)
		}
		val tailLen = args(0).toInt
		val mRnaFile = args(1)
		val outputFile = args(2)
		val basesPerLine = if (args.length > 3) args(3).toInt else -1
		val polyATail = String.valueOf(Array.fill(tailLen)('A'))
 		val isoformSequences = new IsoformSequencesFromFile(mRnaFile, null, false)
		println("writing isoform sequences to " + outputFile + "...");
		val pw = new PrintWriter(outputFile);
		for (isoform <- asIterable(isoformSequences.getAllTags)) {
			val isoSequence = isoformSequences.getSequence(isoform) + polyATail
			pw.write('>')
			pw.write(isoform.toString)
			pw.write('\n')

			if (basesPerLine == -1) {
				pw.write(isoSequence);
				pw.write('\n');
			} else {
				for (i <- 0 until(isoSequence.length, basesPerLine)){
					pw.write(isoSequence.substring(i, min(i + basesPerLine, isoSequence.length)))
					pw.write('\n');
				}
			}
		}
		pw.close();
	}
}


