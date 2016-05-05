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
package edu.uconn.engr.dna.util;

import java.util.Arrays;

/**
 *
 * @author marius
 */
public class TrieJ {

	private TrieNode root;
	private int[] index;
	private char[] alphabet;

	public TrieJ(char[] alphabet) {
		this.alphabet = alphabet;
		this.root = new TrieNode(new TrieNode[alphabet.length]);
		this.index = createChildIndex(alphabet);
	}

	public void add(String path, String isoform, int pos) {
		TrieNode node = root;
		int n = path.length() - 1;
		for (int i = 0; i < n; ++i) {
			node = node.childM(index[path.charAt(i)]);
		}
		node.leafChildM(index[path.charAt(n)]).add(isoform, pos);
	}

	class Finder {

		private String path;
		PrimitiveStringArray strings;
		PrimitiveIntArray ints;
		boolean copied = false;

		public Finder(String path) {
			this.path = path;
		}

		public Pair<String[], int[]> find(int mismatches) {
			find(root, 0, mismatches);
			if (strings != null) {
				return new Pair<String[], int[]>(strings.toArray(), ints.toArray());
			} else {
				return null;
			}
		}

		private void find(TrieNode root, int pathIndex, int mismatches) {
			if (root == null) {
				return;
			} else if (pathIndex == path.length()) {
				TrieLeaf leaf = (TrieLeaf) root;
				if (strings == null) {
					strings = leaf.strings;
					ints = leaf.ints;
				} else {
					if (!copied) {
						strings = strings.plus(leaf.strings);
						ints = ints.plus(leaf.ints);
						copied = true;
					} else {
						strings.addAll(leaf.strings);
						ints.addAll(leaf.ints);
					}
				}
			} else {
				char c = path.charAt(pathIndex);
				int k = index[c];
				if (k >= 0) {
					// first, exact match
					find(root.child(k), pathIndex + 1, mismatches);
				}
				if (mismatches > 0) {
					for (int i = 0; i < alphabet.length; ++i) {
						if (i != k) {
							TrieNode child = root.child(i);
							if (child != null) {
								// then, allow mismatches
								find(child, pathIndex + 1, mismatches - 1);
							}
						}
					}
				}
			}
		}
	}

	public Pair<String[], int[]> find(final String path, int mismatches) {
		return new Finder(path).find(mismatches);
	}

	private int[] createChildIndex(char[] alphabet) {
		int[] ind = new int[256];
		Arrays.fill(ind, -1);
		for (int i = 0; i < alphabet.length; ++i) {
			ind[alphabet[i]] = i;
		}
		return ind;
	}

	/*	public int totalNodes() {
	return root.totalNodes();
	}

	public int totalLeaves() {
	return root.totalLeaves();
	}

	public int totalElements() {
	return root.totalElements();
	} */
	class TrieNode {

		private TrieNode[] children;

		public TrieNode(TrieNode[] children) {
			this.children = children;
		}

		public TrieNode childM(int index) {
			if (children[index] == null) {
				return children[index] = new TrieNode(new TrieNode[children.length]);
			}
			return children[index];
		}

		public TrieNode child(int index) {
			return children[index];
		}

		public TrieLeaf leafChildM(int index) {
			if (children[index] == null) {
				return (TrieLeaf) (children[index] = new TrieLeaf());
			}
			return (TrieLeaf) children[index];
		}

		public boolean isLeaf() {
			return false;
		}
		/*
		public int nElements {
		return 0;
		}

		def totalNodes: Int = if (isLeaf) 1 else children.foldLeft(1){(s, c) =>
		if (c == null) s else s + c.totalNodes
		}

		def totalLeaves: Int = if (isLeaf) 1 else children.foldLeft(0){(s, c) =>
		if (c == null) s else s + c.totalLeaves
		}

		def totalElements: Int = if (isLeaf) nElements else children.foldLeft(0){(s, c) =>
		if (c == null) s else s + c.totalElements
		} */
	}

	class TrieLeaf extends TrieNode {

		private PrimitiveStringArray strings = new PrimitiveStringArray();
		private PrimitiveIntArray ints = new PrimitiveIntArray();

		public TrieLeaf() {
			super(null);
		}

		public void add(String s, int i) {
			strings.add(s);
			ints.add(i);
		}

		@Override
		public boolean isLeaf() {
			return true;
		}
	}
}
