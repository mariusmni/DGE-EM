package edu.uconn.engr.dna.dgeem

import java.util.Arrays.asList
import joptsimple.{OptionSet, OptionParser}

/**
 * User: marius
 * Date: Oct 12, 2010
 */

class DGEOptionParser extends OptionParser {
  import DGEOptionParser._
  acceptsAll(asList(OP_HELP, "help"), "Show help")
  acceptsAll(asList(OP_GTF, "GTF"), "Known genes and isoforms in GTF format")
  		.withRequiredArg.describedAs("GTF file")
  acceptsAll(asList(OP_CLUSTERS_FILE, "gene-clusters"), 
		  "Override isoform to gene mapping defined in the " 
		+ "GTF file with a mapping taken from the given file. The format " 
		+ "of each line in the file is \"isoform\tgene\".")
		.withRequiredArg.describedAs("cluster file")
  accepts(OP_MISMATCHES, "Maximum number of mismatched allowed for a tag (default 1)." ).withRequiredArg.ofType(classOf[java.lang.Integer])
//  acceptsAll(asList(OP_QUALITY, "quality-scores"),
//		  "Weigh the tags based on their quality scores.")
  accepts(OP_UNIQ, "Infer frequencies only from tags that map to the same gene")
//  accepts(OP_PROB, "Cleave probability").withRequiredArg.ofType(classOf[java.lang.Double])
  accepts(OP_MRNA, "Isoform MRNA sequences in 5' to 3' orientation").withRequiredArg
  accepts(OP_ENZIME, "Enzime cutting patterns (comma separated, no spaces)").withRequiredArg.describedAs("enzymes")
  accepts(OP_PREPEND_ENZYME, "Use this flag if the recognition site for the enzyme is not included in the tags")
  accepts(OP_OUTPUT_PREFIX, "Output files prefix (default 'dge')")
  		.withRequiredArg.describedAs("prefix")
  accepts(OP_LIMIT_NTAGS, "Discard all tags after this many have been read").withRequiredArg.ofType(classOf[java.lang.Integer])
}

object DGEOptionParser {
//  final val OP_QUALITY = "q"
  final val OP_HELP = "h"
  final val OP_MISMATCHES = "max-mismatches"
  final val OP_CLUSTERS_FILE = "c"
  final val OP_GTF = "G"
  final val OP_MRNA = "m"
  final val OP_UNIQ = "uniq"
//  final val OP_PROB = "p"
  final val OP_ENZIME = "e"
  final val OP_PREPEND_ENZYME = "prepend-enzyme"
  final val OP_LIMIT_NTAGS = "limit-ntags"
  final val OP_OUTPUT_PREFIX = "o"

  def main(args: Array[String]) {
	  new DGEOptionParser().printHelpOn(System.out)
  }
}