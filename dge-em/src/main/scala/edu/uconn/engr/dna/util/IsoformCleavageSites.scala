
package edu.uconn.engr.dna.util;
import edu.uconn.engr.dna.format.Clusters
import edu.uconn.engr.dna.format.Isoforms

import scala.collection.{mutable => cm}
import scala.collection.mutable.ListBuffer;
import scala.collection.JavaConversions._;

import edu.uconn.engr.dna.format.TaggedSequences

object IsoformCleavageSites {
	val patternMatch = Array.fill[Array[Char]](256)(null)
	for (p <- List(
			'R' -> Array('G', 'A'),
			'Y' -> Array('C', 'T'),
			'M' -> Array('A', 'C'),
			'K' -> Array('G', 'T'),
			'S' -> Array('G', 'C'),
			'W' -> Array('A', 'T'),
			'B' -> Array('C', 'G', 'T'),
			'D' -> Array('A', 'G', 'T'),
			'H' -> Array('A', 'C', 'T'),
			'V' -> Array('A', 'C', 'G'),
			'N' -> Array('A', 'C', 'G', 'T'))) {
		patternMatch(p._1) = p._2
	}


	def println(m: String) {
		Console.err.println("sites:" + m)
	}

	def alphabet(ts: TaggedSequences) = {
		val set = cm.Set[Char]()
		for (tag <- ts.getAllTags;
			 c <- ts.getSequence(tag).toString) {
			set.add(c)
		}
		set.toArray
	}

	def checkIsoforms(ts: TaggedSequences, isoforms: Isoforms, clusters: Clusters): Boolean = {
		println("checking isoforms consistence...")
		ts.getAllTags.filterNot(t => isoforms.contains(t.toString)).foreach(ts.remove(_))

		val toRemove = for{
			id <- isoforms.idIterator;
			sequence = ts.getSequence(id) if (sequence == null)
				} yield id

			if (toRemove.size > 0) {
				println("Erroneous isoforms " + toRemove.size)
				val isoformsToClustersMap = Utils.createIsoformToClusterMap(isoforms, clusters)
				toRemove.foreach(id => {
						isoforms.remove(id)
						val cluster = clusters.get(isoformsToClustersMap.get(id))
						cluster.remove(id)
						if (cluster.isEmpty)
							clusters.remove(cluster.getName)
						println("ERROR: no sequence could be found for isoform " + id + "! isoform removed")
					})
			}
			toRemove.size > 0
		}

		/**
		 * Same as the two parameter version, except that it discards positions
		 * closer than tagLen to the end of the isoform
		 */
		def sitesForIsoforms(isoformSequences: TaggedSequences, cleavePatterns: List[String], tagLen: Int) =
			computeIsoformCleaveSites(isoformSequences, cleavePatterns, tagLen).map(p =>
				(p._1, p._2 .size));
	
		def sitesForIsoforms(isoformSequences: TaggedSequences, cleavePatterns: List[String]) =
			computeIsoformCleaveSites(isoformSequences, cleavePatterns).map(p => (p._1, p._2 .size));

		/**
		 * Same as the two parameter version, except that it discards positions
		 * closer than tagLen to the end of the isoform
		 */
		def computeIsoformCleaveSites(isoformSequences: TaggedSequences,
									  cleavePatterns: List[String], tagLen: Int): Map[String, IndexedSeq[Int]] = {
			computeIsoformCleaveSites(isoformSequences, cleavePatterns).map(p =>{
					val isoLength = isoformSequences.getSequence(p._1).length
					(p._1, p._2.filter(pos => pos + tagLen <= isoLength))
				})
		}

		def computeIsoformCleaveSites(isoformSequences: TaggedSequences,
									  cleavePatterns: List[String]): Map[String, IndexedSeq[Int]] = {
			val map = for{
				isoformName <- isoformSequences.getAllTags;
				s = isoformSequences.getSequence(isoformName).toString
			} yield (isoformName.toString, cleavePatterns.toIndexedSeq[String].
					 flatMap(pattern => getAllMatchPositions(pattern, s)).
					 sortWith(_ > _))
			map.toMap
		}

		def getAllMatchPositions(pattern: String, s: String): Iterable[Int] = {
			val positions = ListBuffer[Int]()
			var pos = indexOf(pattern, s, 0)
			while (pos >= 0) {
				positions.add(pos)
				pos = indexOf(pattern, s, pos + 1)
			}
			positions
		}

		def indexOf(pattern: String, s: String, pos: Int): Int = {
			val limit = (s.length - pattern.length)
			var i = pos
			while (i < limit) {
				var j=0
				while (j < pattern.length && matches(pattern.charAt(j), s.charAt(i+j))) {
					j+=1
				}
				if (j == pattern.length) {
					return i
				}
				i += 1
			}
			return -1
		}

		def matches(patternSymbol: Char, sequenceSymbol: Char) = {
			(patternSymbol == sequenceSymbol) ||
			(patternMatch(patternSymbol) != null &&
			 patternMatch(patternSymbol).contains(sequenceSymbol))
		}
	}