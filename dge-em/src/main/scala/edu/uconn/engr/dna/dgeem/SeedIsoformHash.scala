package edu.uconn.engr.dna.dgeem;
import edu.uconn.engr.dna.util.IsoformCleavageSites;
import edu.uconn.engr.dna.format.TaggedSequences;
import edu.uconn.engr.dna.util.Utils
import scala.collection.JavaConversions._;
import scala.collection.{mutable => cm};

/**
 * Class for finding matches between tags and isoforms and computing
 * weights based on quality scores (with mismatches).
 */
class SeedIsoformHash(isoformSequences: TaggedSequences, cleavePatterns: List[String]) {
	var fresh = true

	var map: Map[String, TagClass] = null

	var count = 0
	var countKeys = 0
	val dummyTag = new TagClass
	
	def getTagClass(tag: String, qual: String) = {
		if (fresh) {
			this.synchronized {
				if (fresh) { 
					map = initialize(tag.length);
					println("Total sites with mismatches " + map.size)
					//	println("Total isoforms matched " + map.foldLeft(0){
					//			case (s, (tag, tc)) => s + tc.length
					//		})
					println("Total isoforms matched " + count)
					fresh = false
				}
			}
		};
		val o = map.get(tag)
		if (o != None) {
			val tc = o.get
			new WeightedTagClass(tc.iso.toArray, tc.pos.toArray, weights(tag, qual, tc).toArray, 1)
		} else {
			null
		}
	}

	def initialize(tagLen: Int) = {
		println("Initializing sites");
		val sitesForIsoform = IsoformCleavageSites
		.computeIsoformCleaveSites(isoformSequences,
								   cleavePatterns, tagLen);
		println("Total sites " + (sitesForIsoform.foldLeft(0)((s, l) => s + l._2.size)))
		
		tagClassesForSites(tagLen, sitesForIsoform, atMostOneMismatch(identity))//atMostOneMismatch(identity)))
	};

	def tagClassesForSites(tagLen: Int,
						   sitesForIsoform: Map[String, IndexedSeq[Int]],
						   alteredTags: (String) => Seq[String] = identity) = {
		val map = collection.mutable.Map[String, TagClass]()
		var isos = 0
		for ((isoform, sites) <- sitesForIsoform;
		     sequence = isoformSequences.getSequence(isoform).toString) {
			isos += 1
			for (i <- 0 until sites.length;
				 site = sites(i) if (site + tagLen <= sequence.length);
				 tag = sequence.substring(site, site + tagLen);
				 alteredTag <- alteredTags(tag)) {
				hash(map, alteredTag)//.addMatch(isoform, i)
				count+=1
				if (countKeys % 100000 == 0) {
					println("Key count " + countKeys + " isos processed " + isos + " total neighbors processed " + count)
				}
			}
		}
		map.toMap
	}
	
	def hash(map: cm.Map[String, TagClass], tag: String) = {
		val o = map.get(tag)
		if (o == None) {
			val nt = dummyTag // new TagClass()
			countKeys += 1
			map(tag) = nt
			nt
		} else {
			o.get
		}
	}

	def atMostOneMismatch(tailTransform: (String) => Seq[String])(tag: String) = {
		val n = new cm.ListBuffer[String]
		n.add(tag)
		for (i <- 12 until tag.length;
			 c = tag.charAt(i) if (c != 'N');
			 before = tag.substring(0, i);
			 after <- tailTransform(tag.substring(i+1))) {
			var d = nextBase(c);
			while (d != c) {
				n.add(before + d + after)
				d = nextBase(d)
			}
		}
		println("Neighbors of " + tag + " " + n.length + ": " + n)
		n
	}

	def identity(tag: String) = List(tag)

	def nextBase(b: Char) = b.toUpper match {
		case 'A' => 'C'
		case 'C' => 'T'
		case 'T' => 'G'
		case 'G' => 'A'
	}

	def weights(tag: String, qual: String, tc: TagClass): Array[Double] = {
		val a = new Array[Double](tc.length)
		for (k <- 0 until tc.length;
			 seq = isoformSequences.getSequence(tc.iso(k)).toString;
			 p = tc.pos(k))
				 a(k) = quality(tag, qual, seq, p)
		a
	}

	def quality(tag: String, qual: String, seq: String, pos: Int) = {
		var w = 1.0
		for (i <- 0 until tag.length) {
			val prob = Utils.phredProbability(Utils.fastqToPhred(qual.charAt(i)))
			w *= (if (tag.charAt(i) == seq.charAt(pos+i)) 1-prob else prob/3)
		}
		w
	}
}