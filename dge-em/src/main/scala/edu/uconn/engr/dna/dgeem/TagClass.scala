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

import scala.collection.mutable.ListBuffer

class TagClass(val iso: ListBuffer[String] = new ListBuffer[String],
			   val pos: ListBuffer[Int] = new ListBuffer[Int],
			   var multiplicity: Int = 0) {
	def length = iso.length;
	def addMatch(m: (String, Int)) {
		iso += m._1;
		pos += m._2
	}
	def addMatches(tc: TagClass) {
		iso ++= tc.iso
		pos ++= tc.pos
	}
}


