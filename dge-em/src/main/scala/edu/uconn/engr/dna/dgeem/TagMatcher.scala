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
import edu.uconn.engr.dna.format.TaggedSequences
import edu.uconn.engr.dna.util.IsoformCleavageSites
import edu.uconn.engr.dna.util.Trie
import edu.uconn.engr.dna.util.TrieJ
import edu.uconn.engr.dna.util.Utils
import edu.uconn.engr.dna.util.F
import collection.{mutable => cm}
import java.{util => ju}
import scala.collection.{JavaConversions => jc}
import scala.collection.mutable.Builder

class TagMatcher(isoformSequences: TaggedSequences, 
				 cleavePatterns: List[String],
				 mismatches: Int,
				 qualScores: Boolean,
				 enzymeIncluded: Boolean) {
	var noMatch = 0
	var fresh = true
	private var sitesForIsoform: ju.Map[String, Array[Int]] = null
	var nrCleaveSites: Map[String, Int] = null
	
	var trie: TrieJ = null

	def getTagClass(tag: String, qual: String) = {
		if (fresh) {
			this.synchronized {
				if (fresh) {
					initialize(tag.length);
					fresh = false
				}
			}
		};
		val p = trie.find(tag, mismatches)
		if (p == null) {
//			println("Could not find a match for " + tag)
			noMatch += 1
			null
		} else {
			val isoArray = p.getFirst
			val posArray = p.getSecond
			val w = if (qualScores) 
				F.normalize(weights(tag, qual, isoArray, posArray))
			else
				Array.fill(isoArray.length)(1.0/isoArray.length)
			new WeightedTagClassJ(isoArray, posArray, w, 1)
		}
	}

	def initialize(tagLen: Int) = {
		println("Initializing sites");
		val sites = IsoformCleavageSites.computeIsoformCleaveSites(isoformSequences,
								   cleavePatterns, tagLen);
		sitesForIsoform = new ju.HashMap[String, Array[Int]]
		if (!enzymeIncluded) {
			val offset = cleavePatterns.head.length
			for ((k,v) <- sites) {
				sitesForIsoform.put(k, v.map(x => x + offset).toArray)
			}
		} else {
			for ((k,v) <- sites) {
				sitesForIsoform.put(k, v.toArray)
			}
		}
		println("Total sites " + (sites.foldLeft(0)((s, l) => s + l._2.size)))
		nrCleaveSites = sites.map(p => (p._1, p._2 .size))

		tagClassesForSites(tagLen, jc.asScalaMap(sitesForIsoform))
	};


	def tagClassesForSites(tagLen: Int,
						   sitesForIsoform: cm.Map[String, Array[Int]]) {
		val alphabet = cm.Set[Char]()
		for ((isoform, sites) <- sitesForIsoform;
		     sequence = isoformSequences.getSequence(isoform).toString;
			 i <- 0 until sites.length;
			 site = sites(i) if (site + tagLen <= sequence.length);
			 tag = sequence.substring(site, site + tagLen);
			 c <- tag) {
			alphabet.add(c)
		}
		println("Alphabet: " + alphabet)
		trie = new TrieJ(alphabet.toArray)
		val sitesSortedByIsoform = sitesForIsoform.toArray.sortWith((p, q) => p._1.compareTo(q._1) <= 0)
		for ((isoform, sites) <- sitesSortedByIsoform;
		     sequence = isoformSequences.getSequence(isoform).toString;
			 i <- 0 until sites.length;
			 site = sites(i) if (site + tagLen <= sequence.length);
			 tag = sequence.substring(site, site + tagLen)) {
			trie.add(tag, isoform, i)
		}
	}

	def weights(tag: String, qual: String, isos: Array[String], siteInds: Array[Int]):Array[Double] = {
		val result = new Array[Double](isos.length)
		var i = 0
		while (i < isos.length) {
			val iso = isos(i)
			val siteIndex = siteInds(i)
			val seq = isoformSequences.getSequence(iso).toString
			val pos = sitesForIsoform.get(iso)(siteIndex)
			result(i) =  quality(tag, qual, seq, pos)
			i+=1
		}
		result
	}

	def quality(tag: String, qual: String, seq: String, pos: Int) = {
		var w = 1.0
		var i = tag.length-1
		while (i >= 0) {
			val prob = Utils.phredProbability(Utils.fastqToPhred(qual.charAt(i)))
			w *= (if (tag.charAt(i) == seq.charAt(pos+i)) 1-prob else prob/3)
			i -= 1
		}
		w
	}
}

