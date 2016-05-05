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

class IsoformMatchesCompositeBuffer extends IsoformMatchesBuffer {
	private val strings = new StringArrayWrapper
	private val ints = new IntArrayWrapper

	override def += (e: (String, Int)) {
		strings += e._1
		ints += e._2
	}

	override def ++=(o: IsoformMatchesBuffer) {
		val (s, i) = o.unfold
		strings ++= s
		ints ++= i
	}

	override def length = strings.length 

	override def isEmpty = strings.isEmpty

	override def unfold: (Array[String], Array[Int]) = {
		strings.trimToSize
		ints.trimToSize
		(strings.data, ints.data)
	}

	
}

class StringArrayWrapper {
	var data = new Array[String](8)
	var size = 0
	def length = size
	def isEmpty = size == 0

	def +=(e: String) {
		ensureCapacity(size + 1)
		data(size) = e
		size += 1
	}

	def ++=(a: Array[String]) {
		ensureCapacity(size + a.length)
		Array.copy(a, 0, data, size, a.length)
		size += a.length
	}

	def ensureCapacity(c: Int) {
		if (data.length < c) {
			var k = data.length
			while (k < c)  k <<= 1
			val dest = new Array[String](c)
			Array.copy(data, 0, dest, 0, data.length)
			data = dest
		}
	}

	def trimToSize {
		if (size < data.length) {
			val newData = new Array[String](size)
			Array.copy(data, 0, newData, 0, size)
			data = newData
		}
	}
}

class IntArrayWrapper {
	var data = new Array[Int](8)
	var size = 0
	def length = size
	def isEmpty = size == 0

	def +=(e: Int) {
		ensureCapacity(size + 1)
		data(size) = e
		size += 1
	}

	def ++=(a: Array[Int]) {
		ensureCapacity(size + a.length)
		Array.copy(a, 0, data, size, a.length)
		size += a.length

	}

	def ensureCapacity(c: Int) {
		if (data.length < c) {
			var k = data.length
			while (k < c)  k <<= 1
			val dest = new Array[Int](c)
			Array.copy(data, 0, dest, 0, data.length)
			data = dest
		}
	}

	def trimToSize {
		if (size < data.length) {
			val newData = new Array[Int](size)
			Array.copy(data, 0, newData, 0, size)
			data = newData
		}
	}
}