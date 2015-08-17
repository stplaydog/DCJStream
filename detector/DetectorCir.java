package detector;

import graphs.Graph;
import tools.Func;

public class DetectorCir extends Detector {

	public void detect_ASs(Graph g, int up_num) {
		this.clean();
		if (this.AS1(g))
			return;
		else if (AS2_one(g)) {
			if (this.num_detected > 1) {
				this.transMajor();
				if (this.AS4(g))
					return;
				else
					this.transMajorBack();
			} else
				return;
		} else if (this.AS4(g))
			return;
		else
			this.AS0(g);
		return;
	}

	private boolean AS1(Graph g) {
		for (int i = 0; i < g.node_num; i++)
			valid[i] = g.check[i];
		for (int v = 0; v < g.node_num; v++) {
			if (!valid[v])
				continue; // when v is incident to a selected edge
			vn = g.adj_mat[v];
			if (vn[0] == vn[1] || vn[0] == vn[2]) {
				this.addVertex(v);
				this.addVertex(vn[0]);
				valid[v] = valid[vn[0]] = false;
			} else if (vn[1] == vn[2]) {
				this.addVertex(v);
				this.addVertex(vn[1]);
				valid[v] = valid[vn[1]] = false;
			}
		}
		// end of search

		if (this.idx_major > 0) { // find AS1
			this.num_detected = 1;
			return true;
		} else
			return false;
	}

	private boolean AS2_one(Graph g) {
		boolean two = false; // whether a four cycle with only two colors has
								// been found
		// here we use two boolean arrays: valid0 and valid
		// valid0 record whether a vertex is incident to a selected edge
		// valid record whether a vertex belongs to an examed subgraph
		for (int i = 0; i < g.node_num; i++)
			valid0[i] = g.check[i];

		// loop1: we check for each color
		for (int color = 0; color < 3; color++) {
			for (int i = 0; i < g.node_num; i++)
				valid[i] = valid0[i];

			// loop2: we check each vertex under a given color
			// 'i' is the vertex we start
			// 'j' is its adj_mat via the edge of the given color
			for (int i = 0; i < g.node_num; i++) {
				if (!valid[i])
					continue; // the vertex is invalid here if it is in a
								// discovered AS subgraph or has been searched
								// already.
				int j = g.adj_mat[i][color];
				if (j >= g.node_num || !valid[j])
					continue;
				valid[i] = false;
				valid[j] = false;

				int c1 = (color + 1) % 3, c2 = (color + 2) % 3; // c1 and c2 are
																// the other two
																// colors

				// the four vertice adj_mating to i and j
				int i1 = g.adj_mat[i][c1], i2 = g.adj_mat[i][c2];
				int j1 = g.adj_mat[j][c1], j2 = g.adj_mat[j][c2];

				if (g.adj_mat[i1][color] == j2 && valid[i1] && valid[j2]) {
					// 0 1' 2, 0 --- each number represents an edge of that
					// color
					this.addVertex(i);
					this.addVertex(i1);
					this.addVertex(j);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i1] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i1] = valid[j2] = false;
				} else if (g.adj_mat[i2][color] == j1 && valid[i2] && valid[j1]) {
					// 0 2' 1, 0
					this.addVertex(i);
					this.addVertex(i2);
					this.addVertex(j);
					this.addVertex(j1);
					valid0[i] = valid0[j] = valid0[i2] = valid0[j1] = false;
					valid[i] = valid[j] = valid[i2] = valid[j1] = false;
				} else if (g.adj_mat[i1][c2] == j1 && valid[i1] && valid[j1]) {
					// 0 1' 1, 2
					this.addVertex(i);
					this.addVertex(j);
					this.addVertex(i1);
					this.addVertex(j1);
					valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
					valid[i] = valid[j] = valid[i1] = valid[j1] = false;
				} else if (g.adj_mat[i2][c1] == j2 && valid[i2] && valid[j2]) {
					// 0 2' 2, 1

					this.addVertex(i);
					this.addVertex(j);
					this.addVertex(i2);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i2] = valid[j2] = false;
				}

				else if (g.adj_mat[i1][color] == j1 && (i2 == j1 || j2 == i1)
						&& valid[i1] && valid[j1]) {
					// 0 1' 1, 0 i2 == j1 || j2 == i1
					this.addVertex(i);
					this.addVertex(i1);
					this.addVertex(j);
					this.addVertex(j1);
					valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
					valid[i] = valid[j] = valid[i1] = valid[j1] = false;
				} else if (g.adj_mat[i2][color] == j2 && (i1 == j2 || j1 == i2)
						&& valid[i2] && valid[j2]) {
					// 0 2' 2, 0 i1 == j2 || j1 == i2
					this.addVertex(i);
					this.addVertex(i2);
					this.addVertex(j);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i2] = valid[j2] = false;
				} else if (g.adj_mat[i1][color] == j1
						&& g.adj_mat[i2][color] == j2 && i1 != j2 && i2 != j1
						&& valid[i1] && valid[i2] && valid[j1] && valid[j2]) {
					// double four 0 1' 1, 0 2' 2, 0
					this.addVertex(i);
					this.addVertex(j);
					this.addVertex(i1);
					this.addVertex(j1);
					this.addVertex(i2);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i1] = valid0[i2] = valid0[j1] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i1] = valid[i2] = valid[j1] = valid[j2] = false;
				} else if (this.idx_major == 0 && two == false) {
					if (g.adj_mat[i1][color] == j1 && valid[i1] && valid[j1]) {
						// System.out.println("it's c!");
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i1;
						four_cycle[3] = j1;
					} else if (g.adj_mat[i2][color] == j2 && valid[i2]
							&& valid[j2]) {
						// System.out.println("it's d!");
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i2;
						four_cycle[3] = j2;
					}
				}
			} // end of loop 2
		} // end of loop 1

		if (this.idx_major > 0 || two) {
			if (this.idx_major > 0) {
				this.num_detected = 1;
				return true;
			} else {
				this.addVertex(four_cycle[0]);
				this.addVertex(four_cycle[1]);
				this.addVertex(four_cycle[2]);
				this.addVertex(four_cycle[3]);
				this.addVertex(four_cycle[0]);
				this.addVertex(four_cycle[2]);
				this.addVertex(four_cycle[1]);
				this.addVertex(four_cycle[3]);
				this.num_detected = 2;
				return true;
			}
		} else
			return false;
	}

	// public boolean AS2_all(Graph g) {
	// Vector<Vector<Integer>> two = new Vector<Vector<Integer>>();
	// Vector<Integer> t = new Vector<Integer>();
	// int four_cycle[] = new int[4]; // to store a four cycle which only
	// // contains two colors
	// // boolean two = false; // whether a
	// // four cycle with only two colors has
	// // been found
	// // here we use two boolean arrays: valid0 and valid
	// // valid0 record whether a vertex is incident to a selected edge
	// // valid record whether a vertex belongs to an examed subgraph
	// boolean valid0[] = new boolean[g.node_num];
	// for (int i = 0; i < g.node_num; i++)
	// valid0[i] = g.check[i];
	//
	// // loop1: we check for each color
	// for (int color = 0; color < 3; color++) {
	// boolean valid[] = new boolean[g.node_num];
	// for (int i = 0; i < g.node_num; i++)
	// valid[i] = valid0[i];
	//
	// // loop2: we check each vertex under a given color
	// // 'i' is the vertex we start
	// // 'j' is its adj_mat via the edge of the given color
	// for (int i = 0; i < g.node_num; i++) {
	// if (!valid[i])
	// continue; // the vertex is invalid here if it is in a
	// // discovered AS subgraph or has been searched
	// // already.
	// int j = g.adj_mat[i][color];
	// if (j >= g.node_num || !valid[j])
	// continue;
	// valid[i] = false;
	// valid[j] = false;
	//
	// int c1 = (color + 1) % 3, c2 = (color + 2) % 3; // c1 and c2 are
	// // the other two
	// // colors
	//
	// // the four vertice adj_mating to i and j
	// int i1 = g.adj_mat[i][c1], i2 = g.adj_mat[i][c2];
	// int j1 = g.adj_mat[j][c1], j2 = g.adj_mat[j][c2];
	//
	// if (g.adj_mat[i1][color] == j2 && valid[i1] && valid[j2]) {
	// // 0 1' 2, 0 --- each number represents an edge of that
	// // color
	//
	// t.add(i);
	// t.add(i1);
	// t.add(j);
	// t.add(j2);
	// valid0[i] = valid0[j] = valid0[i1] = valid0[j2] = false;
	// valid[i] = valid[j] = valid[i1] = valid[j2] = false;
	// } else if (g.adj_mat[i2][color] == j1 && valid[i2] && valid[j1]) {
	// // 0 2' 1, 0
	// t.add(i);
	// t.add(i2);
	// t.add(j);
	// t.add(j1);
	// valid0[i] = valid0[j] = valid0[i2] = valid0[j1] = false;
	// valid[i] = valid[j] = valid[i2] = valid[j1] = false;
	// } else if (g.adj_mat[i1][c2] == j1 && valid[i1] && valid[j1]) {
	// // 0 1' 1, 2
	// t.add(i);
	// t.add(j);
	// t.add(i1);
	// t.add(j1);
	// valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
	// valid[i] = valid[j] = valid[i1] = valid[j1] = false;
	// } else if (g.adj_mat[i2][c1] == j2 && valid[i2] && valid[j2]) {
	// // 0 2' 2, 1
	// t.add(i);
	// t.add(j);
	// t.add(i2);
	// t.add(j2);
	// valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
	// valid[i] = valid[j] = valid[i2] = valid[j2] = false;
	// }
	//
	// else if (g.adj_mat[i1][color] == j1 && (i2 == j1 || j2 == i1)
	// && valid[i1] && valid[j1]) {
	// // 0 1' 1, 0 i2 == j1 || j2 == i1
	// t.add(i);
	// t.add(i1);
	// t.add(j);
	// t.add(j1);
	// valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
	// valid[i] = valid[j] = valid[i1] = valid[j1] = false;
	// } else if (g.adj_mat[i2][color] == j2 && (i1 == j2 || j1 == i2)
	// && valid[i2] && valid[j2]) {
	// // 0 2' 2, 0 i1 == j2 || j1 == i2
	// t.add(i);
	// t.add(i2);
	// t.add(j);
	// t.add(j2);
	// valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
	// valid[i] = valid[j] = valid[i2] = valid[j2] = false;
	// } else if (g.adj_mat[i1][color] == j1
	// && g.adj_mat[i2][color] == j2 && i1 != j2 && i2 != j1
	// && valid[i1] && valid[i2] && valid[j1] && valid[j2]) {
	// // double four 0 1' 1, 0 2' 2, 0
	// t.add(i);
	// t.add(j);
	// t.add(i1);
	// t.add(j1);
	// t.add(i2);
	// t.add(j2);
	// valid0[i] = valid0[j] = valid0[i1] = valid0[i2] = valid0[j1] = valid0[j2]
	// = false;
	// valid[i] = valid[j] = valid[i1] = valid[i2] = valid[j1] = valid[j2] =
	// false;
	// } else if (result.size() == 0) {
	// if (g.adj_mat[i1][color] == j1 && valid[i1] && valid[j1]) {
	// // System.out.println("it's c!");
	// this.addVertex(new Vector<Integer>());
	// result.get(result.size() - 1).add(i); // 0
	// result.get(result.size() - 1).add(j); // 1
	// result.get(result.size() - 1).add(i1); // 2
	// result.get(result.size() - 1).add(j1); // 3
	// this.addVertex(new Vector<Integer>());
	// result.get(result.size() - 1).add(i); // 0
	// result.get(result.size() - 1).add(i1); // 2
	// result.get(result.size() - 1).add(j); // 1
	// result.get(result.size() - 1).add(j1); // 3
	// } else if (g.adj_mat[i2][color] == j2 && valid[i2]
	// && valid[j2]) {
	// // System.out.println("it's d!");
	// this.addVertex(new Vector<Integer>());
	// result.get(result.size() - 1).add(i); // 0
	// result.get(result.size() - 1).add(j); // 1
	// result.get(result.size() - 1).add(i2); // 2
	// result.get(result.size() - 1).add(j2); // 3
	// this.addVertex(new Vector<Integer>());
	// result.get(result.size() - 1).add(i); // 0
	// result.get(result.size() - 1).add(i2); // 2
	// result.get(result.size() - 1).add(j); // 1
	// result.get(result.size() - 1).add(j2); // 3
	// }
	// }
	// } // end of loop 2
	// } // end of loop 1
	// if (result.size() != 0)
	// for (Vector<Integer> res : result)
	// for (int r : t)
	// res.add(r);
	// else {
	// this.addVertex(new Vector<Integer>());
	// for (int r : t)
	// this.addVertex(r);
	// }
	//
	// return true;
	// }

	private boolean AS4(Graph g) {
		for (int i = 0; i < g.node_num; i++)
			valid0[i] = g.check[i];

		/************** step 1: to detect AS4 of type 5-3-5 **************/
		for (int i = 0; i < g.node_num; i++)
			valid[i] = valid0[i];

		out_core_535: for (int core = 0; core < g.node_num; core++) {
			if (!valid[core])
				continue;
			for (int c1 = 0; c1 < 2; c1++) {
				int p1 = g.adj_mat[core][c1];
				if (!valid[p1])
					continue;
				for (int c2 = c1 + 1; c2 < 3; c2++) {
					int p2 = g.adj_mat[core][c2];
					if (!valid[p2])
						continue;
					int c3 = 3 - c1 - c2;
					if (g.adj_mat[p1][c3] == p2) {
						/************* one triangle detected ****************/
						tria[c3] = core;
						tria[c2] = p1;
						tria[c1] = p2;
						valid[core] = valid[p1] = valid[p2] = false;
						// A triangle can be detected three times from its three
						// vertices, but we just need to find it once
						// make these valid false to prevent it being discovered
						// later

						// the three vertices connected to the three endpoints
						// of the triangle; Point Out
						for (int i = 0; i < 3; i++) {
							// to check whether these three vertices are still
							// valid
							po[i] = g.adj_mat[tria[i]][i];
							if (!valid[po[i]])
								continue out_core_535;
						}

						// at this moment we have discovered a valid triangle
						// and its three hands
						/******************* to detect 5-3-5 first ****************/
						// this include subgraphs 2-3,2-4,2-5,2-6
						for (int co1 = 0; co1 < 3; co1++) { // colors extended
															// from Poing Outs
							// we start from the hand with co1 among these three
							// hands
							int co2 = (co1 + 1) % 3, co3 = (co1 + 2) % 3;
							int po12 = g.adj_mat[po[co1]][co2], po13 = g.adj_mat[po[co1]][co3];
							if (!valid[po12] || !valid[po13])
								continue;
							if (g.is_connected(po12, po[co2])
									&& g.is_connected(po13, po[co3])) {
								/************** find a 5-3-5 AS4 *******************/
								this.addVertex(po[co2]);
								this.addVertex(po12);
								this.addVertex(po[co3]);
								this.addVertex(po13);
								this.addVertex(tria[co1]);
								this.addVertex(po[co1]);
								this.addVertex(tria[co2]);
								this.addVertex(tria[co3]);
								valid0[po[co2]] = valid[po[co2]] = false;
								valid0[po12] = valid[po12] = false;
								valid0[po[co3]] = valid[po[co3]] = false;
								valid0[po13] = valid[po13] = false;
								valid0[tria[co1]] = valid[tria[co1]] = false;
								valid0[po[co1]] = valid[po[co1]] = false;
								valid0[tria[co2]] = valid[tria[co2]] = false;
								valid0[tria[co3]] = valid[tria[co3]] = false;
								continue out_core_535; // find a good one and
														// then process next
														// core vertex
							}
						} // end of color
						continue out_core_535;
						// find a triangle but no AS4, process next core vertex;
						// once a vertex is in a triangle, it can not belong to
						// another triangle
					} // end of if it is a triangle
				} // c2
			} // c1
		} // end out_core_535

		/************** step 2: to detect AS4 of type 3-3-3 or 3-3 other types **************/
		for (int i = 0; i < g.node_num; i++)
			valid[i] = valid0[i];

		out_core_333: for (int core = 0; core < g.node_num; core++) {
			if (!valid[core])
				continue; // core may have only two valid adj_mats

			for (int cI = 0; cI < 3; cI++) {
				pI[cI] = g.adj_mat[core][cI];
				for (int cII = 0; cII < 3; cII++) {
					if (!valid[pI[cI]] || cII == cI
							|| !valid[g.adj_mat[pI[cI]][cII]]
							|| !valid[g.adj_mat[g.adj_mat[pI[cI]][cII]][cI]]) {
						pII[cI][cII] = -1;
						pIII[cI][cII] = -1;
					} else {
						pII[cI][cII] = g.adj_mat[pI[cI]][cII];
						pIII[cI][cII] = g.adj_mat[pII[cI][cII]][cI];
					}
				}
			}
			/************* check for 3-3-3 ************/
			for (int cII0 = 0; cII0 < 3; cII0++) {
				if (pIII[0][cII0] == -1)
					continue;
				for (int cII1 = 0; cII1 < 3; cII1++) {
					if (pIII[1][cII1] == -1)
						continue;
					for (int cII2 = 0; cII2 < 3; cII2++) {
						if (pIII[2][cII2] == -1)
							continue;
						if (pIII[0][cII0] == pIII[1][cII1]
								&& pIII[0][cII0] == pIII[2][cII2]) {
							boolean found = false;
							int co_core = pIII[0][cII0];
							/*************** find a cadidate of 3-3-3 **************/
							if (cII0 != cII1 && cII0 != cII2 && cII1 != cII2) {
								found = true;
							}
							/************** find a 3-3-3 ***************/
							else {
								if (g.is_connected(pI[0], pI[1])
										|| g.is_connected(pI[0], pI[2])
										|| g.is_connected(pI[2], pI[1])
										|| g.is_connected(pII[0][cII0],
												pII[1][cII1])
										|| g.is_connected(pII[0][cII0],
												pII[2][cII2])
										|| g.is_connected(pII[2][cII2],
												pII[1][cII1])) {
									found = true;
									/************ find a 3-3x3 ***************/
								}
							}
							if (found) {

								this.addVertex(core);
								this.addVertex(co_core);
								this.addVertex(pI[0]);
								this.addVertex(pII[0][cII0]);
								this.addVertex(pI[1]);
								this.addVertex(pII[1][cII1]);
								this.addVertex(pI[2]);
								this.addVertex(pII[2][cII2]);
								valid0[core] = valid[core] = false;
								valid0[co_core] = valid[co_core] = false;
								valid0[pI[0]] = valid[pI[0]] = false;
								valid0[pII[0][cII0]] = valid[pII[0][cII0]] = false;
								valid0[pI[1]] = valid[pI[1]] = false;
								valid0[pII[1][cII1]] = valid[pII[1][cII1]] = false;
								valid0[pI[2]] = valid[pI[2]] = false;
								valid0[pII[2][cII2]] = valid[pII[2][cII2]] = false;
								continue out_core_333;
							} // end of found
						} // end of 333
					} // end of cII2
				} // end of cII1
			} // end of cII0

			// at this stage, we can find at most 3-3 pattern
			/************* check for 3-3-other ************/
			for (int c1 = 0; c1 < 2; c1++) { // start from core with color c1
				for (int c2 = c1 + 1; c2 < 3; c2++) { // start from core with
														// color c2
					for (int c1c = 0; c1c < 3; c1c++) { // start from pI[c1]
														// with color c1c
						if (pIII[c1][c1c] == -1)
							continue;
						for (int c2c = 0; c2c < 3; c2c++) { // start from pI[c2]
															// with c2c
							if (pIII[c2][c2c] == -1)
								continue;
							if (pIII[c1][c1c] == pIII[c2][c2c]) {
								/******** find a 3-3 **********/
								int co_core = pIII[c1][c1c];
								int out1 = -1, out2 = -1;
								boolean has_33_other = false;
								int c3 = 3 - c1 - c2;
								if (valid[g.adj_mat[core][c3]]
										&& valid[g.adj_mat[co_core][c3]]) {
									// if there are valid "three hands"-- the
									// vertices incident to core/co_core via the
									// third color
									out1 = g.adj_mat[core][c3];
									out2 = g.adj_mat[co_core][c3];
									if (out1 == g.adj_mat[g.adj_mat[core][c1]][3
											- c1 - c1c]
											&& out2 == g.adj_mat[g.adj_mat[co_core][c1]][3
													- c1 - c1c]) {
										has_33_other = true;
										// find an instance of 3-3 or 3-4
									} else if (out1 == g.adj_mat[g.adj_mat[core][c2]][3
											- c2 - c2c]
											&& out2 == g.adj_mat[g.adj_mat[co_core][c2]][3
													- c2 - c2c]) {
										has_33_other = true;
										// find an instance of 3-3 or 3-4
									}
								}
								if (!has_33_other && c1 == c2c && c2 == c1c) { // a
																				// cycle
																				// of
																				// size
																				// 6
									/************ try to detect (3-5) **************/
									int p1 = g.adj_mat[core][c1], p2 = g.adj_mat[core][c2];
									int cop1 = g.adj_mat[co_core][c1], cop2 = g.adj_mat[co_core][c2];
									// p1, p2, cop1, cop2 are guranteed to be
									// valid
									if (g.adj_mat[p1][c3] == p2) {
										out1 = g.adj_mat[cop1][c3];
										out2 = g.adj_mat[cop2][c3];
										if (valid[out1] && valid[out2]
												&& g.is_connected(out1, out2)) {
											has_33_other = true;
											// find an instance of (3-5)
										}
									} else if (g.adj_mat[cop1][c3] == cop2) {
										out1 = g.adj_mat[p1][c3];
										out2 = g.adj_mat[p2][c3];
										if (valid[out1] && valid[out2]
												&& g.is_connected(out1, out2)) {
											has_33_other = true;
											// find an instance of (3-5)
										}
									}
								}
								if (!has_33_other) {
									/*********** try to detect (3-1) or (3-2) **************/
									int p1 = g.adj_mat[core][c1], p2 = g.adj_mat[core][c2];
									int cop1 = g.adj_mat[co_core][c1], cop2 = g.adj_mat[co_core][c2];
									int p1e = g.adj_mat[p1][3 - c1 - c1c], cop2e = g.adj_mat[cop2][3
											- c2 - c2c];
									int p2e = g.adj_mat[p2][3 - c2 - c2c], cop1e = g.adj_mat[cop1][3
											- c1 - c1c];
									out1 = p1e;
									out2 = p2e;
									if (valid[out1] && out1 == cop2e
											&& valid[out2] && out2 == cop1e) {
										has_33_other = true;
										// find an instance of (3-1) or (3-2)
									}
								}

								// final to record
								if (has_33_other) {
									this.addVertex(core);
									this.addVertex(co_core);
									this.addVertex(out1);
									this.addVertex(out2);
									this.addVertex(pI[c1]);
									this.addVertex(pII[c1][c1c]);
									this.addVertex(pI[c2]);
									this.addVertex(pII[c2][c2c]);
									valid0[core] = valid[core] = false;
									valid0[co_core] = valid[co_core] = false;
									valid0[out1] = valid[out1] = false;
									valid0[out2] = valid[out2] = false;
									valid0[pI[c1]] = valid[pI[c1]] = false;
									valid0[pII[c1][c1c]] = valid[pII[c1][c1c]] = false;
									valid0[pI[c2]] = valid[pI[c2]] = false;
									valid0[pII[c2][c2c]] = valid[pII[c2][c2c]] = false;
									continue out_core_333;
								}
							} // find a 3-3 structure
						} // c2c
					} // c1c
				} // c2
			} // c1
		} // out_core_333
		if (this.idx_major == 0)
			return false;
		this.num_detected = 1;
		return true;
	}

	private void AS0(Graph g) {

		int point = 0;
		if (!this.is_zero) {
			// use the heuristics to minimize the search space or maximize the
			// possible lower bound for heuristics
			g.get_start(0, 1);
			g.get_start(0, 2);
			g.get_start(1, 2);

			if (this.is_heu)
				point = g.cal_low_start();
			else
				point = g.cal_up_start();
		} else {
			// just select the point with lowest id
			for (int i = 0; i < g.node_num; i++)
				if (g.check[i]) {
					point = i;
					break;
				}
		}
		// add all possible combinations
		for (int i = 0; i < g.node_num; i++) {
			if (g.check[i] != false && i != point) {
				this.addVertex(point);
				this.addVertex(i);
				this.num_detected++;
			}
		}
	}

}
