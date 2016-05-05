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
import edu.uconn.engr.dna.io.GTFParser
import java.io.FileInputStream

object TestTagMatcher {
	def main(args: Array[String]) {
		println("Parsing GTF file...");
		val (clusters, isoforms) = parseGTF("knownGeneGnfAtlas2.gtf")
		println("Found " + isoforms.size + " isoforms and " + clusters.size + " genes");

		println("reading known gene sequences...");
		// Isoforms sequences are loaded as they are in the file,
		// which is in coding strand orientation!
		val isoformSequencesFile = "myKnownGeneMrna.txt";
		val isoformSequences = new IsoformSequencesFromFile(isoformSequencesFile, isoforms, false);

		val tm = new TagMatcher(isoformSequences, List("CATG"), mismatches=2, qualScores=true, enzymeIncluded=true)
		println(tm.getTagClass("ACTGACTGACTGACTGACTGA", "[[[[[[[[[[[[[[[[[[[[["))
	}

	def parseGTF(file: String): Tuple2[Clusters, Isoforms] = {
		val giParser = new GTFParser;
		val p = giParser.parse(new FileInputStream(file));
		(p.getFirst, p.getSecond)
	}
}
