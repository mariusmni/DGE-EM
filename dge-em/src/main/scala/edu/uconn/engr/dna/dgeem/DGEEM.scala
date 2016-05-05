package edu.uconn.engr.dna.dgeem
import collection.mutable.Set
import edu.cornell.lassp.houle.RngPack.Ranlux
import edu.uconn.engr.dna.format.Isoforms
import collection.breakOut
import java.{util => ju, lang => jl}
import collection.{JavaConversions => jc}

/**
 * User: marius
 * Date: Sep 2, 2010
 */

class DGEEM(var p: Double, isoforms: Isoforms, nrCleavageSites: ju.Map[String, Int]) {

	private val maxSteps = 500
	private val EPS = 0.000000000001
	private val EPS2 = EPS * EPS

	private val random = new Ranlux(4, 123);

	import DGEEM._

	def computeFreq(tagClasses: Array[WeightedTagClassJ]): ju.Map[String, Double] = {
		p = 0.1
		println("Initialize p to " + p)

		val isoformsSet = new ju.HashSet[String]
		for (tc <- tagClasses) {
			for (i <- tc.isoforms) {
				isoformsSet.add(i)
			}
		}
		val uniqIsoforms = new Array[String](isoformsSet.size)
		var i = 0
		for (isoform <- jc.asScalaIterable(isoformsSet)) {
			uniqIsoforms(i) = isoform
			i+=1
		}
		val n = uniqIsoforms.size
		
		val indexJ = new ju.HashMap[jl.String, jl.Integer]
		var currentIndex = 0
		for (i <- uniqIsoforms) {
			indexJ.put(i, currentIndex)
			currentIndex += 1
		}

		val nrCleaveSites = new Array[Int](n)
		uniqIsoforms.foreach(i => nrCleaveSites(indexJ.get(i).intValue) = nrCleavageSites.get(i))

		val count = new Array[Array[Double]](n)
		for (i <- 0 until n) count(i) = new Array[Double](nrCleaveSites(i))
		
		val maxSites = nrCleaveSites.max 

		var f = Array.fill(uniqIsoforms.size)(1.0)
//		for (i <- uniqIsoforms) f.put(i, 1.0 /*random.uniform(0.0001, 1.0)*/)
		var fSum = f.sum

		// speed up future operations by converting isoform names to indices
		i = tagClasses.length-1
		while (i >= 0) {
			val tc = tagClasses(i)
			tc.indexData(indexJ)
			i-=1
		}

		var steps = 0
		var maxDiff = 1.0
		while (maxDiff > 0.00005 && steps < maxSteps) {
			val q = 1 - p
			val qp = preprocessQPowers(q, maxSites)

			// E step
			var i = tagClasses.length-1
			while (i >= 0) {
				val tc = tagClasses(i)
				i-=1
				var sum = 0.0

				var k = tc.size-1
				while (k >= 0) {
					val i = tc.iso(k)
					val j = tc.pos(k)
					val w = tc.w(k)
					k -= 1
					sum +=  f(i) * qp(j) * w //! make sure j ranges from 0 to nrCleavage(i)-1
				}

				if (sum > EPS2) {
					k = tc.size-1
					while (k >= 0) {
						val i = tc.iso(k)
						val j = tc.pos(k)
						val w = tc.w(k)
						k -= 1
						val fraction = f(i) * qp(j) * w / sum
						val increment = tc.multiplicity * fraction
//						println(i + " gets fraction " + fraction + " of " +
//										tagClass.multiplicity + " = " + increment)
						count(i)(j) += increment
					}
				} else {
//					println("Give up on " + tc.multiplicity + " tags because sum is " + sum)
				}
			}

			// M step
			val newFreq = new Array[Double](n)
			for (i <- 0 until n) {
				var r = (1-qp(nrCleaveSites(i)))
				if (r < EPS) r = 1.0;
				newFreq(i) = count(i).sum / r
			}

			// re-estimate p
			var N1 = 0.0
			var D = 0.0
			for (array <- count) {
				N1 += array(0)
				D += array.sum
			}
			p = N1/D
			println("Number of tags aligning to first position " + N1)
			println("Number of tags aligning to any position " + D)
			println("p reestimated to " + p)


			// convergence test
			val newFSum = newFreq.sum
			maxDiff = 0.0
			for (i <- 0 until n) {
				val oldNormVal = f(i) / fSum;
				if (oldNormVal > EPS) {
					val relDiff = ((newFreq(i) / newFSum) - oldNormVal) / oldNormVal
					if (relDiff > maxDiff) maxDiff = relDiff
				}
			}

			// prepare for next iteration
			f = newFreq
			fSum = newFSum
			fillWithZeros(count)
			steps += 1
//			println("Step " + steps + " classes " + tagClasses.size)
		}
		if (steps == maxSteps) {
			lock.synchronized {
				noConverge+=1
			}
		}
		println("p converged to " + p)
		val freq = new ju.HashMap[String, Double]
		for (i <- 0 until n) {
			freq.put(uniqIsoforms(i), f(i))
		}
		freq
	}

	def findMax[E](items: Iterable[E], map: ju.Map[E, Int]) = {
		var m = 0
		for (i <- items) {
			if (map.containsKey(i)) {
				val v = map.get(i)
				if (v > m) {
					m = v
				}
			}
		}
		m
	}
	
	def preprocessQPowers(q: Double, maxSites: Int) = {
		val qp = new Array[Double](maxSites+1)
		qp(0) = 1
		for (i <- 1 to maxSites) qp(i) = qp(i-1) * q
		qp
	}
  
	def power(a: Double, b: Int): Double = {
		if (b == 0)
			1
		else if ((b & 1) == 0)
			power(a * a, b >> 1)
		else
			a * power(a * a, b >> 1)
	}

	def fillWithZeros(m: Array[Array[Double]]) {
		for (a <- m) {
			fillWithZeros(a)
		}
	}
	
	def fillWithZeros(a: Array[Double]) = for (i <- 0 until a.size) a(i) = 0.0

	def sumValues(m: ju.Map[String, Double]) = {
		val iterator = m.entrySet.iterator
		var sum = 0.0
		while (iterator.hasNext) {
			sum += iterator.next.getValue
		}
		sum
	}

}

object DGEEM {
	val lock = new Object()
	var noConverge = 0
}