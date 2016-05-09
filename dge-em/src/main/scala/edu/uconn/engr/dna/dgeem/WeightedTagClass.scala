/*
 *  Copyright 2010 marius.
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


import java.{util => ju}
import scala.collection.{JavaConversions => jc}
import java.lang.{Double => JDouble}

class WeightedTagClass(var isoforms: Array[String],
					   val pos: Array[Int],
					   val w: Array[Double],
					   var multiplicity: Int = 0) {
	assert(isoforms.length == pos.length && pos.length == w.length, 
		   "All arrays must have the same length; found " + isoforms.length +
		   " " + pos.length + " " + w.length)

	var iso : Array[Int] = null
	def size = isoforms.length
	
	def incrementMultiplicity() { multiplicity += 1}
	
	def indexData(index: Map[String, Int]) {
		iso = isoforms.map(i => index(i))
 	}

	def normalizeWeights {
		val sum = w.sum
		for (i <- 0 until w.length) {
			w(i) /= sum
		}
	}

	def sort {
		val sorted = (isoforms, pos, w).zipped.toList.sortWith((a, b) => {
					val c = a._1.compareTo(b._1) // compare isoform names
					c < 0 || (c == 0 && a._2 < b._2) // if equal, compare positions
				})
		var k = 0
		for ((i, p, v) <- sorted) {
			isoforms(k) = i
			pos(k) = p
			w(k) = v
			k += 1
		}
	}

	override def equals(o: Any): Boolean = {
		if (!o.isInstanceOf[WeightedTagClass]) {
			return false
		}
		val t = o.asInstanceOf[WeightedTagClass]
        if (size != t.size) {
            return false
		}

        for (i <- 0 until isoforms.length) {
            if (pos(i) != t.pos(i) ||
				JDouble.doubleToLongBits(w(i)) != JDouble.doubleToLongBits(t.w(i)) ||
				isoforms(i) != t.isoforms(i)) {
                return false
			}
		}
		return true
	}

	override def hashCode = {
		var h = 0
		for (i <- 0 until isoforms.length) {
			h = 31 * h + pos(i)
            val bits = JDouble.doubleToLongBits(w(i))
			h = 31 * h + (bits ^ (bits >>> 32)).asInstanceOf[Int]
			val element = isoforms(i)
			h = 31 * h + (if (element == null) 0 else element.hashCode())
		}
		h
	}
}


