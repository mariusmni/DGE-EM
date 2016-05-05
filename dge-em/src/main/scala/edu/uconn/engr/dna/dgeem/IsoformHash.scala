package edu.uconn.engr.dna.dgeem;
import edu.uconn.engr.dna.util.IsoformCleavageSites;
import edu.uconn.engr.dna.format.TaggedSequences;
import scala.collection.JavaConversions._;
import scala.collection.{mutable => cm};

class IsoformHash(isoformSequences: TaggedSequences, cleavePatterns: List[String]) {
	var fresh = true

	var map: Map[String, WeightedTagClass] = null
	var sitesForIsoform: Map[String, IndexedSeq[Int]] = null

	def getAllNonZeroTagClasses = map.values.filter(_.multiplicity > 0)
	
	def incrementMatches(tag: String) {
		if (fresh) {
			this.synchronized {
				if (fresh) { 
					val tmap = initialize(tag.length);
					map = tmap.map{case (k, v) =>
							(k -> new WeightedTagClass(v.iso.toArray,
													   v.pos.toArray,
													   Array.fill(v.iso.size)(1.0)))}

					fresh = false
				}
			}
		};
		val o = map.get(tag)
		if (o != None) {
			val tc = o.get
			tc.synchronized { 
				tc.multiplicity+=1
			}
		}
	}

	def initialize(tagLen: Int) = {
		val map = collection.mutable.Map[String, TagClass]()
		println("Initializing sites");
		sitesForIsoform = IsoformCleavageSites.computeIsoformCleaveSites(isoformSequences, 
										     cleavePatterns, tagLen);
		for ((isoform, sites) <- sitesForIsoform;
		     sequence = isoformSequences.getSequence(isoform).toString;
		     i <- 0 until sites.length;
		     site = sites(i) if (site + tagLen <= sequence.length);
		     tag = sequence.substring(site, site + tagLen)) {
			hash(map, tag, (isoform, i))
		}
		println("Total sites " + (sitesForIsoform.foldLeft(0)((s, l) => s + l._2.size)))
		map.toMap
	};
	
	def nrCleaveSites = sitesForIsoform.map(p => (p._1, p._2 .size))

	def hash(map: cm.Map[String, TagClass], tag: String, t: (String, Int)) {
		val o = map.get(tag)
		val tc = if (o == None) {
			val nt = new TagClass()
			map(tag) = nt
			nt
		} else {
			o.get
		}
		tc.addMatch(t)
	}
}