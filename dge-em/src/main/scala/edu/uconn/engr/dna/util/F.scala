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

import scala.collection.generic.CanBuildFrom

object F {
	def zip[A, B, C, That](a: Iterable[A], b: Iterable[B], fn: (A,B)=>C)(implicit bf: CanBuildFrom[Iterable[A], C, That]): That = {
		val builder = bf(a)
		val these = a.iterator
		val those = b.iterator
		while (these.hasNext && those.hasNext)
			builder += fn(these.next, those.next)
		builder.result
	}

	def zip3[A, B, C, D, That](a: Iterable[A], b: Iterable[B], c: Iterable[C], fn: (A,B,C)=>D)
		(implicit bf: CanBuildFrom[Iterable[A], D, That]): That = {
		val builder = bf(a)
		val these = a.iterator
		val those = b.iterator
		val those2 = c.iterator
		while (these.hasNext && those.hasNext && those2.hasNext)
			builder += fn(these.next, those.next, those2.next)
		builder.result
	}

	def factory[A,B](thunk: () => ParameterRunnable[A, B]) = new ParameterRunnableFactory[A,B]() {
    	override def createParameterRunnable = thunk()
	}

	def converter[A,B](conv: A=>B) = new Converter[A,B]() {
		override def convert(a: A) = conv(a)
	}

	def binaryOperator[A, B, C](fn: (A,B)=>C) = new BinaryOperator[A, B, C]() {
		override def compute(a: A, b: B) = fn(a, b)
	}

	def normalize(a: Array[Double]): Array[Double] = {
		var s = 0.0
		var i = 0
		while (i < a.length) {
			s += a(i)
			i += 1
		}
		i = 0
		while (i < a.length) {
			a(i) /= s
			i += 1
		}
		a
	}
}
