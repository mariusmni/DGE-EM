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

package edu.uconn.engr.dna.util
import java.{util => ju, lang => jl}

class CompositeBuffer[E] extends IsoformMatchesBuffer {
	val components: ju.List[IsoformMatchesBuffer] = new ju.ArrayList[IsoformMatchesBuffer]()

	override def += (e: (String, Int)) {
		throw new UnsupportedOperationException("Operation not supported")
	}

	override def ++= (o: IsoformMatchesBuffer) {
		components.add(o)
	}

	override def unfold: (Array[String], Array[Int]) = {
		val n = length
		val nstr = new Array[String](n)
		val nint = new Array[Int](n)
//		var i = n
		var i = 0
		var j = components.size-1
		while (j >= 0) {
			val bb = components.get(j)
			j -= 1
			val (bbs, bbi) = bb.unfold
//			merge(bbs, bbi, 0,
//				  nstr, nint, i,
//				  nstr, nint, i-bbs.length)
			System.arraycopy(bbs, 0, nstr, i, bbs.length)
			System.arraycopy(bbi, 0, nint, i, bbi.length)
//			i -= bbs.length
			i += bbs.length
		}
		(nstr, nint)
	}
	
	def merge(as: Array[String], ai: Array[Int], af: Int,
			  bs: Array[String], bi: Array[Int], bf: Int,
			  os: Array[String], oi: Array[Int], of: Int) {
		var i = af
		var j = bf
		var k = of
		while (i < as.length && j < bs.length) {
			var c = as(i).compareTo(bs(j))
			if (c == 0) {
				c = ai(i) - bi(j)
			}
			if (c < 0) {
				os(k) = as(i)
				oi(k) = ai(i)
				i+=1
			} else {
				os(k) = bs(j)
				oi(k) = bi(j)
				j+=1
			}
			k += 1
		}
		while (i < as.length) {
			os(k) = as(i)
			oi(k) = ai(i)
			i+=1
			k += 1
		}
		while (j < bs.length) {
			os(k) = bs(j)
			oi(k) = bi(j)
			j+=1
			k += 1
		}
	}
		
	override def length: Int = if (components == null) 0 else {
		var s = 0
		var i = components.size - 1
		while (i >= 0) {
			s += components.get(i).length
			i -= 1
		}
		s
	}

	override def isEmpty = components.isEmpty
}
