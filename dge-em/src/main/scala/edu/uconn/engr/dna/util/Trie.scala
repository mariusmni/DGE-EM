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

import java.{util => ju}
import scala.collection.{JavaConversions => JC}
import scala.collection.{mutable => cm}


class Trie(alphabet: Array[Char]) {
	type E = (String, Int)
	val root = new TrieNode(new Array[TrieNode](alphabet.size));
	val index = createChildIndex(alphabet)
	val weirdos = cm.Set[Char]()
	alphabet.foreach(e => weirdos.add(e))
	
	def add(path: String, element: E) {
		var node = root
		for (i <- 0 until (path.length-1)) {
			node = node.child(index(path.charAt(i)), true)
		}
		node.leafChild(index(path.last), true).add(element)
	}

	def find(path: String, mismatches: Int): IsoformMatchesBuffer = {
		def find(root: TrieNode, pathIndex: Int, mismatches: Int): IsoformMatchesBuffer = {
			if (root == null) {
				null
			} else if (pathIndex == path.length) {
				root.asInstanceOf[TrieLeaf].buffer
			} else {
				val c = path.charAt(pathIndex)
				val k = index(c)
				var result: IsoformMatchesBuffer = if (k >= 0) {
					// first, exact match
					find(root.child(k, false), pathIndex+1, mismatches)
				} else {
					this.synchronized {
						if (!weirdos.contains(c)) {
							weirdos.add(c)
							println("Unexpected character " + c)
						}
					}
					null
				}
				if (mismatches > 0) {
					var copied = false
					var i = alphabet.length-1
					while (i >= 0) {
						if (i != k) {
							val child = root.child(i, false)
							if (child != null) {
								// then, allow mismatches
								val m = find(child, pathIndex+1, mismatches-1)
								if (m != null) {
									if (result == null) {
										result = m
									} else  {
										if (!copied) {
											val newBuffer = new IsoformMatchesCompositeBuffer
											newBuffer ++= result
											result = newBuffer
											copied = true
										}
										result ++= m
									}
								}
							}
						}
						i -= 1
					}
				}
				result
			}
		}
		find(root, 0, mismatches)
	}

	private def createChildIndex(alphabet: Array[Char]) = {
		val index = Array.fill(256)(-1)
		for (i <- 0 until alphabet.length) {
			index(alphabet(i)) = i
		}
		index
	}


	def totalNodes = root.totalNodes

	def totalLeaves = root.totalLeaves

	def totalElements = root.totalElements


	class TrieNode(children: Array[TrieNode]) {
		def child(index: Int, create: Boolean) = {
			if (create && children(index) == null) {
				children(index) = new TrieNode(new Array[TrieNode](children.length))
			}
			children(index)
		}

		def leafChild(index: Int, create: Boolean): TrieLeaf = {
			if (create && children(index) == null) {
				children(index) = new TrieLeaf(new IsoformMatchesCompositeBuffer)
			}
			children(index).asInstanceOf[TrieLeaf]
		}

		def isLeaf = false
		def nElements = 0

		def totalNodes: Int = if (isLeaf) 1 else children.foldLeft(1){(s, c) =>
			if (c == null) s else s + c.totalNodes
		}

		def totalLeaves: Int = if (isLeaf) 1 else children.foldLeft(0){(s, c) =>
			if (c == null) s else s + c.totalLeaves
		}

		def totalElements: Int = if (isLeaf) nElements else children.foldLeft(0){(s, c) =>
			if (c == null) s else s + c.totalElements
		}
	}

	class TrieLeaf(val buffer: IsoformMatchesBuffer) extends TrieNode(null) {
		def add(e: E) {
			buffer += e
		}
		override def isLeaf = true
		override def nElements = buffer.length
	}
}
