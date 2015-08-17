package gasts;


import java.util.ArrayList;

public class WeightedASM {
	double cap_thr1=1; // the threshold for AS1 when the cap is involved
//		double cap_thr2=0.875; // for large edge length >=0.875 edge=400
//	double cap_thr2=1; // for large edge length >=0.875 edge=400
		double cap_thr2=0.85; // for small edge length >=0.6 edge=380
	int num_gene;
	int num_chr;
	int size;
	int num_free_caps;
	Edges e1, e2;
	WeightedAdjacency wa;
	Adjacency median;
	ArrayList<Integer> available;

	public WeightedASM (AdjNode n1, AdjNode n2, AdjNode n3, AdjNode avoid) {		
		if(n1.adj==null || n2.adj==null) {
			System.out.println("In WeightedASM the first two nodes are not initialized!");
			System.exit(1);
		}

		e1=new Edges(1,n1); 
		e2=new Edges(1,n2);
		wa=new WeightedAdjacency();

		ArrayList<AdjNode> searched_already=new ArrayList<AdjNode>();
		searched_already.add(avoid);
		ArrayList<AdjNode> to_search=new ArrayList<AdjNode>();
		ArrayList<AdjNode> next_to_search=new ArrayList<AdjNode>();
		to_search.add(n3);

		int dist=0;
		int max_num_chr=Math.max(n1.adj.num_chr, n2.adj.num_chr);

		while(!to_search.isEmpty()) {
			float weight=(float) Math.pow(2, -dist);
			for(AdjNode current_node:to_search) {
				searched_already.add(current_node);

				if(current_node.isLeaf) {
					wa.add_leaf(weight, current_node);
					if(current_node.adj.num_chr>max_num_chr) {
//						System.out.println(current_node.name+" "+current_node.adj.num_chr);
						max_num_chr=current_node.adj.num_chr;
					}
				}

				for(int i=0;i<current_node.links.size();i++) {
					AdjNode n=(AdjNode) current_node.links.get(i);
					if(n==null) continue;
					if(!searched_already.contains(n)) next_to_search.add(n);
				}
			}
			to_search=next_to_search;
			next_to_search=new ArrayList<AdjNode>();
			dist++;
		}

		median=new Adjacency(n1.adj.num_gene,max_num_chr, "median");
		num_gene=median.num_gene;
		num_chr=median.num_chr;
		size=num_gene+num_chr;
		num_free_caps=2*num_chr;
		available=new ArrayList<Integer>(2*num_gene);
		for(int i=0;i<2*num_gene;i++) available.add(i);

		wa.update();
		//		wa.print();
	}

	public void get_median() {
		ArrayList<Integer> result;		
		do{
			//			check();
			result=as1();
			//			System.out.println("as1 "+size);
			if(result==null || result.isEmpty()) {
				result=as2();
				//				System.out.println("as2 "+size);
			}
			//			if(result==null || result.isEmpty()) {
			//				result=as2more();
			//				System.out.println("as2more "+size);
			//			}
			if(result==null || result.isEmpty()) {
				result=within_cycles(true);
//								System.out.println("within_cycles "+size);
			}
//			if(result==null || result.isEmpty()) {
//				result=as2_2wa();
//				if(!result.isEmpty()) System.out.println("as2_2wa "+size+" "+result.size());
////				if(result.size()>0) System.out.println("as2_2wa "+size);
//			}
//			if(result==null || result.isEmpty()) {
//				result=conserved();
//				//				System.out.println("conserved "+size);
//			}
//			if(result==null || result.isEmpty()) {
//				result=within_cycles(false);
//								System.out.println("within_cycles any"+size);
//			}

			if(result==null || result.isEmpty()) {
				result=ash();
				//				System.out.println("ash "+size);
			}
			if(!result.isEmpty()) add_median_edges(result);
//			else {
//				result=new ArrayList<Integer>(available);
//				if(available.size()%2==1) result.add(Constant.CAP0);
//				add_median_edges(result);
//			}
			//						e1.print();
			//						e2.print();
			//									wa.print();
			//									System.out.println();
//			if(available.size()+num_free_caps!=2*size) 
//				System.out.println();
		}while(available.size()>0 && !result.isEmpty());


//		if(size>0) 
//			System.out.println("Hasn't finished yet! "+size );

		median.num_chr=median.cap_adj.size()/2;
//				 System.out.println();
	}

	ArrayList<Integer> as1 () {
		int cap_cnt=num_free_caps;
		ArrayList<Integer> result=new ArrayList<Integer>();
		boolean valid[]=new boolean[available.size()];
		for(int i=0;i<available.size();i++) valid[i]=true;

		// on e1
		for(int left:available) {
			int right=e1.regular[left];
			if(right<left || right>=Constant.NULL) continue;
			if(right==e2.regular[left]) {
				if(right<Constant.CAP1) {
					result.add(left);result.add(right);
					valid[index(left)]=false;
					valid[index(right)]=false;
				}
				else if(right==Constant.CAP0 && cap_cnt>0) {
					cap_cnt--;
					result.add(left);result.add(right);
					valid[index(left)]=false;
				}
			}
			else if(right<Constant.CAP1 && wa.regular[left][right]>=0.5) {	
				result.add(left);result.add(right);
				valid[index(left)]=false;
				valid[index(right)]=false;
			}
			else if(right==Constant.CAP0 && wa.cap0[left]>=cap_thr1 && cap_cnt>0) {
				cap_cnt--;
				result.add(left);result.add(right);
				valid[index(left)]=false;
			}
		}

		// on e2 
		for(int left:available) {
			if(!valid[index(left)]) continue;
			int right=e2.regular[left];
			if(right<left || right>=Constant.NULL) continue;
			if(right<Constant.CAP1 && wa.regular[left][right]>=0.5) {	
				if(!valid[index(right)]) continue;
				result.add(left);result.add(right);
				valid[index(left)]=false;
				valid[index(right)]=false;
			}
			else if(right==Constant.CAP0 && wa.cap0[left]>=cap_thr1 && cap_cnt>0) {
				cap_cnt--;
				result.add(left);result.add(right);
				valid[index(left)]=false;
			}
		}

		// cap AS
		for(int l:e1.cap0) {
			if(!valid[index(l)]) continue;
			for(int r:e1.cap1) {
				if(!valid[index(r)]) continue;
				if(wa.regular[l][r]>=cap_thr2) {
					valid[index(l)]=valid[index(r)]=false;
					result.add(l);result.add(r);
				}
			}
		}
		return result;
	}

	//	ArrayList<Integer> as1p () {
	//		// to find multiple edges between e1,wa or e2,wa
	//		int cap_cnt=num_free_caps;
	//		ArrayList<Integer> result=new ArrayList<Integer>();
	//		boolean[] valid=new boolean[available.size()];
	//		for(int i=0;i<available.size();i++) valid[i]=true;
	//
	//		// on e1
	//		for(int left: available) {
	//			int i_l=index(left);
	//			if(!valid[i_l]) continue;
	//			int right=e1.regular[left];
	//			if(right>=Constant.NULL) continue;
	//			if(right<Constant.CAP1) {
	//				int i_r=index(right);
	//				if(!valid[i_r]) continue;
	//				if(wa.regular[left][right]>=0.5) {
	//					valid[i_l]=valid[i_r]=false;
	//					result.add(left);result.add(right);
	//				}
	//				else if(wa.cap0[left]+wa.cap1[right]>=1 || wa.cap1[left]+wa.cap0[right]>=1) {
	//					valid[i_l]=valid[i_r]=false;
	//					result.add(left);result.add(right);
	//				}
	//			}
	//			else if(right==Constant.CAP0 && wa.cap0[left]>=0.5 && cap_cnt>0) {
	//				valid[i_l]=false;
	//				result.add(left);result.add(Constant.CAP0);
	//				cap_cnt--;
	//			}
	//		}
	//
	//		// cap AS
	//		for(int l:e1.cap0) 
	//			for(int r:e1.cap1) {
	//				if(valid[index(l)] && valid[index(r)] && wa.regular[l][r]>=0.5) {
	//					valid[index(l)]=valid[index(r)]=false;
	//					result.add(l);result.add(r);
	//				}
	//			}
	//
	//		// on e2
	//		for(int left:available) {
	//			int i_l=index(left);
	//			if(!valid[i_l]) continue;
	//			int right=e2.regular[left];
	//			if(right>=Constant.NULL) continue;
	//			if(right<Constant.CAP1) {
	//				int i_r=index(right);
	//				if(!valid[i_r]) continue;
	//				if(wa.regular[left][right]>=0.5) {
	//					valid[i_l]=valid[i_r]=false;
	//					result.add(left);result.add(right);
	//				}
	//				else if(wa.cap0[left]+wa.cap1[right]>=1 || wa.cap1[left]+wa.cap0[right]>=1) {
	//					valid[i_l]=valid[i_r]=false;
	//					result.add(left);result.add(right);
	//				}
	//			}
	//			else if(right==Constant.CAP0 & wa.cap0[left]>=0.5 &&  cap_cnt>0) {
	//				valid[i_l]=false;
	//				result.add(left);result.add(Constant.CAP0);
	//				cap_cnt--;
	//			}
	//
	//		}
	//		for(int l:e2.cap0) 
	//			for(int r:e2.cap1) {
	//				if(valid[index(l)] && valid[index(r)] && wa.regular[l][r]>=0.5) {
	//					valid[index(l)]=valid[index(r)]=false;
	//					result.add(l);result.add(r);
	//				}
	//			}
	//		return result;
	//	}

	ArrayList<Integer> as2p () {
		// to find multiple edges between e1,wa or e2,wa
		ArrayList<Integer> result=new ArrayList<Integer>();
		boolean[] valid=new boolean[available.size()];
		for(int i=0;i<available.size();i++) valid[i]=true;

		// on e1
		for(int left: available) {
			int i_l=index(left);
			if(!valid[i_l]) continue;
			int right=e1.regular[left];
			if(right>=Constant.CAP1) continue;
			int i_r=index(right);
			if(!valid[i_r]) continue;
			int left2=e2.regular[left], right2=e2.regular[right];
			if(left2>=Constant.CAP1 || right2>=Constant.CAP1) continue;
			int i_l2=index(left2), i_r2=index(right2);
			if(!valid[i_l2] || !valid[i_r2]) continue;
			if(wa.regular[left2][right2]==1) {
				valid[i_l]=valid[i_r]=false;
				valid[i_l2]=valid[i_r2]=false;
				result.add(left);result.add(right);
				result.add(left2);result.add(right2);
			}
		}

		// on e1
		for(int left: available) {
			int i_l=index(left);
			if(!valid[i_l]) continue;
			int right=e2.regular[left];
			if(right>=Constant.CAP1) continue;
			int i_r=index(right);
			if(!valid[i_r]) continue;
			int left2=e1.regular[left], right2=e1.regular[right];
			if(left2>=Constant.CAP1 || right2>=Constant.CAP1) continue;
			int i_l2=index(left2), i_r2=index(right2);
			if(!valid[i_l2] || !valid[i_r2]) continue;
			if(wa.regular[left2][right2]==1) {
				valid[i_l]=valid[i_r]=false;
				valid[i_l2]=valid[i_r2]=false;
				result.add(left);result.add(right);
				result.add(left2);result.add(right2);
			}
		}
		return result;
	}

	ArrayList<Integer> conserved () {
		ArrayList<Integer> result=new ArrayList<Integer>();
		boolean[] valid=new boolean[available.size()];
		for(int i=0;i<available.size();i++)  valid[i]=true;

		// on wa
		for(int i=0;i<available.size();i++) {
			if(!valid[i]) continue;
			for(int j=i+1;j<available.size();j++) {
				if(!valid[j]) continue;
				int left=available.get(i), right=available.get(j);
				if(wa.regular[left][right]==1) {
					valid[i]=valid[j]=false;
					result.add(left);result.add(right);
				}
			}
		}

		if(num_free_caps>0) {
			int caps_cnt=num_free_caps;
			for(int i=0;i<available.size();i++) {
				if(!valid[i]) continue;
				int left=available.get(i);
				if(wa.cap0[left]==1 && caps_cnt>0) {
					valid[i]=false;
					result.add(left);result.add(Constant.CAP0);
					caps_cnt--;
				}
			}
		}

		return result;
	}

	ArrayList<Integer> as1h () {

		float max_w=-1;
		int the_left=0, the_right=0;
		int cap_cnt=num_free_caps;
		ArrayList<Integer> result=new ArrayList<Integer>();
		boolean[] valid=new boolean[available.size()];
		for(int i=0;i<available.size();i++)  valid[i]=true;

		// on wa
		for(int i=0;i<available.size();i++) {
			if(!valid[i]) continue;
			for(int j=i+1;j<available.size();j++) {
				if(!valid[j]) continue;
				int left=available.get(i), right=available.get(j);
				if(wa.regular[left][right]==1) {
					valid[i]=valid[j]=false;
					result.add(left);result.add(right);
				}
			}
		}

		if(num_free_caps>0) {
			int caps_cnt=num_free_caps;
			for(int i=0;i<available.size();i++) {
				if(!valid[i]) continue;
				int left=available.get(i);
				if(wa.cap0[left]==1 && caps_cnt>0) {
					valid[i]=false;
					result.add(left);result.add(Constant.CAP0);
					caps_cnt--;
				}
			}
		}

		if(result.size()>0) 
			return result;

		// on e1
		for(int left:available) {
			int right=e1.regular[left];
			if(right>=Constant.NULL || right<left) continue;
			if(right<Constant.CAP1) {
				if(wa.regular[left][right]>max_w) {
					max_w=wa.regular[left][right];
					the_left=left;
					the_right=right;
				}
			}
			else if(right==Constant.CAP0 && cap_cnt>0 && wa.cap0[left]>max_w) {
				max_w=wa.cap0[left];
				the_left=left;
				the_right=right;
			}
		}

		// on e2
		for(int left=0;left<2*num_gene;left++) {
			int right=e2.regular[left];
			if(right>=Constant.NULL || right<left) continue;
			if(right<Constant.CAP1) {
				if(wa.regular[left][right]>max_w) {
					max_w=wa.regular[left][right];
					the_left=left;
					the_right=right;
				}
			}
			else if(right==Constant.CAP0 && cap_cnt>0 && wa.cap0[left]>max_w) {
				max_w=wa.cap0[left];
				the_left=left;
				the_right=right;
			}
		}

		//		System.out.println("max weight="+max_w+" size is "+size);
		if(max_w>=0) {
			result.add(the_left);result.add(the_right);
		}
		else if(max_w<0) result=new ArrayList<Integer>(available);
		return result;
	}

	ArrayList<Integer> ash () {	
		float max_w=-1;
		int cap_cnt=num_free_caps;
		ArrayList<Integer> result=new ArrayList<Integer>();

		// on e1
		for(int left:available) {
			int right=e1.regular[left];
			if(right>=Constant.NULL || right<left) continue;
			if(right<Constant.CAP1) {
				float w=(float) (wa.regular[left][right]-0.5);
				if(w>max_w) {
					max_w=w;
					result=new ArrayList<Integer>();
					result.add(left);
					result.add(right);
				}
			}
			else if(right==Constant.CAP0 && cap_cnt>0) { 
				float w=(float) (wa.cap0[left]-0.5);
				if(w>max_w) {
					max_w=w;
					result=new ArrayList<Integer>();
					result.add(left);
					result.add(right);
				}
			}
		}

		// on e2
		for(int left=0;left<2*num_gene;left++) {
			int right=e2.regular[left];
			if(right>=Constant.NULL || right<left) continue;
			if(right<Constant.CAP1) {
				float w=(float) (wa.regular[left][right]-0.5);
				if(w>max_w) {
					max_w=w;
					result=new ArrayList<Integer>();
					result.add(left);
					result.add(right);
				}
			}
			else if(right==Constant.CAP0 && cap_cnt>0) { 
				float w=(float) (wa.cap0[left]-0.5);
				if(w>max_w) {
					max_w=w;
					result=new ArrayList<Integer>();
					result.add(left);
					result.add(right);
				}
			}
		}


		// cycles of 4
		// e1
		for(int left: available) {
			//			if(!valid[i_l]) continue;
			int right=e1.regular[left];
			if(right>=Constant.CAP1) continue;
			//				if(!valid[i_r]) continue;
			int left2=e2.regular[left], right2=e2.regular[right];
			if(left2>=Constant.CAP1 || right2>=Constant.CAP1) continue;
			//				if(!valid[i_l2] || !valid[i_r2]) continue;
			float w=(wa.regular[left2][right2]-1)/2;
			if(w>max_w) {
				//					valid[i_l]=valid[i_r]=false;
				//					valid[i_l2]=valid[i_r2]=false;
				max_w=w;
				result=new ArrayList<Integer>();
				result.add(left);result.add(right);
				result.add(left2);result.add(right2);
			}
		}

		// on e1
		for(int left: available) {
			//			if(!valid[i_l]) continue;
			int right=e2.regular[left];
			if(right>=Constant.CAP1) continue;
			//			if(!valid[i_r]) continue;
			int left2=e1.regular[left], right2=e1.regular[right];
			if(left2>=Constant.CAP1 || right2>=Constant.CAP1) continue;
			//			if(!valid[i_l2] || !valid[i_r2]) continue;
			float w=(wa.regular[left2][right2]-1)/2;
			if(w>max_w) {
				max_w=w;
				//				valid[i_l]=valid[i_r]=false;
				//				valid[i_l2]=valid[i_r2]=false;
				result=new ArrayList<Integer>();
				result.add(left);result.add(right);
				result.add(left2);result.add(right2);
			}
		}

		//				System.out.println("max weight="+max_w+" size is "+size+" size of the chocie "+result.size());

		if(max_w>=-0.5) {
			//			System.out.println(result.get(result.size()-1)+" "+result.get(result.size()-2));
			return result;
		}
		else if(max_w<-0.5) {
			result=new ArrayList<Integer>(available);
			if(available.size()%2==1) result.add(Constant.CAP0);

		}
		return result;
	}

	private int index(int v) { return available.indexOf((Integer) v);}

	ArrayList<Integer> as2 () {
		ArrayList<Integer> result=new ArrayList<Integer>();

		boolean[] valid=new boolean[available.size()]; 
		for(int i=0;i<available.size();i++) valid[i]=true;
		// if only the edges are taken, then the involved vertices are marked as not valid

		// on e1
		for(int left:available) {
			int i_l=index(left);
			if(!valid[i_l]) continue;
			valid[i_l]=false;
			int right=e1.regular[left];
			if(right<Constant.CAP1) {
				int i_r=index(right);
				if(!valid[i_r]) continue;
				valid[i_r]=false;
				int left2=e2.regular[left];
				int right2=e2.regular[right];

				if(left2>=Constant.NULL  || right2>=Constant.NULL ) continue;

				if(left2<Constant.CAP1 && right2<Constant.CAP1){
					int i_l2=index(left2), i_r2=index(right2);
					if(!valid[i_l2] ||!valid[i_r2]) continue;
					if(e1.regular[left2]==right2) {
						// form a 4-cycle
						valid[i_l2]=valid[i_r2]=false;
						float w1=wa.regular[left][right]+wa.regular[left2][right2];
						float w2=wa.regular[left][left2]+wa.regular[right][right2];
						if(w1>=w2) {
							result.add(left);result.add(right);result.add(left2);result.add(right2);
						}
						else {
							result.add(left);result.add(left2);result.add(right);result.add(right2);
						}
					}
					else if(wa.regular[left][right]+wa.regular[left2][right2]>=1.0) {
						valid[i_l2]=valid[i_r2]=false;
						result.add(left);result.add(right);result.add(left2);result.add(right2);
					}
				}
				else if(left2==Constant.CAP0 && right2==Constant.CAP1 || left2==Constant.CAP1 && right2==Constant.CAP0 ) {
					// a potential 4-cycle with caps involved
					result.add(left);result.add(right);
				}
			}
		}

		// on e2
		for(int left:available) {
			int i_l=index(left);
			if(!valid[i_l]) continue;
			valid[i_l]=false;
			int right=e2.regular[left];
			if(right<Constant.CAP1) {
				int i_r=index(right);
				if(!valid[i_r]) continue;
				valid[i_r]=false;
				int left2=e1.regular[left];
				int right2=e1.regular[right];

				if(left2>=Constant.NULL  || right2>=Constant.NULL ) continue;

				if(left2<Constant.CAP1 && right2<Constant.CAP1){
					int i_l2=index(left2), i_r2=index(right2);
					if(!valid[i_l2] ||!valid[i_r2]) continue;
					if(wa.regular[left][right]+wa.regular[left2][right2]>=1.0) {
						valid[i_l2]=valid[i_r2]=false;
						result.add(left);result.add(right);result.add(left2);result.add(right2);
					}
				}
				else if(left2==Constant.CAP0 && right2==Constant.CAP1 || left2==Constant.CAP1 && right2==Constant.CAP0 ) {
					// a potential 4-cycle with caps involved
					result.add(left);result.add(right);
				}
			}
		}
		return result;
	}

	ArrayList<Integer> as2_2wa() {
		// find a AS2 with 2 edges from wa
		ArrayList<Integer> result=new ArrayList<Integer>();

		// one edge from e1, one edge from e2
		for(int e1_left:available) {
			int e1_right=e1.regular[e1_left];
			if(e1_right<e1_left || e1_right>=Constant.CAP1) continue;
			//			System.out.println();
			for(int e2_left:available) {
				if(e2_left==e1_left || e2_left==e1_right) continue;
				int e2_right=e2.regular[e2_left];
				if(e2_right<e2_left || e2_right>=Constant.CAP1 || e2_right==e1_left || e2_right==e1_right) continue;
				//				System.out.println();
				if(wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_left],wa.regular[e1_right][e2_right])>=1 ||
						wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_right],wa.regular[e1_right][e2_left])>=1) {
					// find a AS2
					//					System.out.println(wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_left],wa.regular[e1_right][e2_right])+" "+
					//							wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_right],wa.regular[e1_right][e2_left]));
					System.out.printf("for e1, e2 edges (%d,%d) %f\t and (%d,%d) %f\n", e1_left, e1_right, wa.regular[e1_left][e1_right], e2_left, e2_right, wa.regular[e2_left][e2_right]);
					System.out.printf("weights\n (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\n\n ", e1_left, e2_left, wa.regular[e1_left][e2_left],
							e1_right, e2_right, wa.regular[e1_right][e2_right],e1_left, e2_right, wa.regular[e1_left][e2_right],e1_right, e2_left, wa.regular[e1_right][e2_left]);
					result.add(e1_left);result.add(e1_right);
					result.add(e2_left);result.add(e2_right);
					System.out.println();
				}

				//				if(wa.regular[e1_left][e2_left]==1 && wa.regular[e1_right][e2_right]==1 || 
				//						wa.regular[e1_left][e2_right]==1 && wa.regular[e1_right][e2_left]==1) {
				//					// find a AS2
				//					result.add(e1_left);result.add(e1_right);
				//					result.add(e2_left);result.add(e2_right);
				//				}
				//				else if(wa.regular[e1_left][e2_left]+wa.regular[e1_right][e2_right]+wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]==2 ||
				//						wa.regular[e1_left][e2_right]+wa.regular[e1_right][e2_left]+wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]==2) {
				//					// find a AS2
				//					result.add(e1_left);result.add(e1_right);
				//					result.add(e2_left);result.add(e2_right);
				//				}
			}
		}

		// two edges from e1
		for(int e1_left:available) {
			int e1_right=e1.regular[e1_left];
			if(e1_right<e1_left || e1_right>=Constant.CAP1) continue;
			//			System.out.println();
			for(int e2_left:available) {
				if(e2_left<=e1_left || e2_left==e1_right) continue;
				int e2_right=e1.regular[e2_left];
				if(e2_right<e2_left || e2_right==e1_right ||e2_right>=Constant.CAP1 ) continue;
				//				System.out.println();
				if(wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_left],wa.regular[e1_right][e2_right])>=1 ||
						wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_right],wa.regular[e1_right][e2_left])>=1) {
					// find a AS2
					System.out.printf("for e1 edges (%d,%d) %f\t and (%d,%d) %f\n", e1_left, e1_right, wa.regular[e1_left][e1_right], e2_left, e2_right, wa.regular[e2_left][e2_right]);
					System.out.printf("weights\n (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\n\n ", e1_left, e2_left, wa.regular[e1_left][e2_left],
							e1_right, e2_right, wa.regular[e1_right][e2_right],e1_left, e2_right, wa.regular[e1_left][e2_right],e1_right, e2_left, wa.regular[e1_right][e2_left]);
					if(wa.regular[e1_left][e2_left]==1 && wa.regular[e1_right][e2_right]==1) {
						result.add(e1_left);result.add(e2_left);
						result.add(e1_right);result.add(e2_right);
					}
					else if(wa.regular[e1_left][e2_right]==1 && wa.regular[e1_right][e2_left]==1) {
						// find a AS2
						result.add(e1_left);result.add(e2_right);
						result.add(e2_left);result.add(e1_right);
					}
					else  {
						// find a AS2
						result.add(e1_left);result.add(e1_right);
						result.add(e2_left);result.add(e2_right);
					}
					System.out.println();
				}

			}
		}

		// two edges from e2

		for(int e1_left:available) {
			int e1_right=e2.regular[e1_left];
			if(e1_right<=e1_left || e1_right>=Constant.CAP1) continue;
			for(int e2_left:available) {
				if(e2_left<=e1_left || e2_left==e1_right) continue;
				int e2_right=e2.regular[e2_left];
				if(e2_right<e2_left || e2_right==e1_right || e2_right>=Constant.CAP1 ) continue;
				if(wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_left],wa.regular[e1_right][e2_right])>=0.9 ||
						wa.regular[e1_left][e1_right]+wa.regular[e2_left][e2_right]+Math.min(wa.regular[e1_left][e2_right],wa.regular[e1_right][e2_left])>=0.9) {
					// find a AS2
					System.out.printf("for e2 edges (%d,%d) %f\t and (%d,%d) %f\n", e1_left, e1_right, wa.regular[e1_left][e1_right], e2_left, e2_right, wa.regular[e2_left][e2_right]);
					System.out.printf("weights\n (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\t (%d,%d): %f\n\n ", e1_left, e2_left, wa.regular[e1_left][e2_left],
							e1_right, e2_right, wa.regular[e1_right][e2_right],e1_left, e2_right, wa.regular[e1_left][e2_right],e1_right, e2_left, wa.regular[e1_right][e2_left]);
					//					result.add(e1_left);result.add(e1_right);
					//					result.add(e2_left);result.add(e2_right);
					if(wa.regular[e1_left][e2_left]==1 && wa.regular[e1_right][e2_right]==1) {
						result.add(e1_left);result.add(e2_left);
						result.add(e1_right);result.add(e2_right);
					}
					else if(wa.regular[e1_left][e2_right]==1 && wa.regular[e1_right][e2_left]==1) {
						// find a AS2
						result.add(e1_left);result.add(e2_right);
						result.add(e2_left);result.add(e1_right);
					}
					else {
						// find a AS2
						result.add(e1_left);result.add(e1_right);
						result.add(e2_left);result.add(e2_right);
					}
					System.out.println();
				}
			}
		}

		return result;
	}

	ArrayList<Integer> as2more () {
		//		boolean remains[]=new boolean[available.size()];
		//		for(int i=0;i<available.size();i++) remains[i]=true;
		//
		//		int e_l, e_r, start;
		//		for(Integer v:e1.cap0) {
		//			if(!remains[index(v)]) continue;
		//			e_r=v;
		//			remains[index(e_r)]=false;
		//
		//			int length=0;
		//			while(true) {
		//				e_l=e2.regular[e_r];
		//				length++;
		//				if(e_l==Constant.CAP0) {System.out.println("00'"+length);break;}
		//				else if(e_l==Constant.CAP1) {System.out.println("01'"+length);break;}
		//				remains[index(e_l)]=false;
		//
		//				e_r=e1.regular[e_l];
		//				length++;
		//				if(e_r==Constant.CAP0) {System.out.println("00"+length);break;}
		//				else if(e_r==Constant.CAP1) {System.out.println("01"+length);break;}
		//				remains[index(e_r)]=false;
		//			}
		//		}
		//
		//		for(Integer v:e1.cap1) {
		//			if(!remains[index(v)]) continue;
		//			e_r=v;
		//			remains[index(e_r)]=false;
		//
		//			int length=0;
		//			while(true) {
		//				e_l=e2.regular[e_r];
		//				length++;
		//				if(e_l==Constant.CAP0) {System.out.println("10'"+length);break;}
		//				else if(e_l==Constant.CAP1) {System.out.println("11'"+length);break;}
		//				remains[index(e_l)]=false;
		//
		//				e_r=e1.regular[e_l];
		//				length++;
		//				if(e_r==Constant.CAP0) {System.out.println("10"+length);break;}
		//				else if(e_r==Constant.CAP1) {System.out.println("11"+length);break;}
		//				remains[index(e_r)]=false;
		//			}
		//		}
		//
		//		for(Integer v:e2.cap0) {
		//			if(!remains[index(v)]) continue;
		//			e_r=v;
		//			remains[index(e_r)]=false;
		//
		//			int length=0;
		//			while(true) {
		//				e_l=e1.regular[e_r];
		//				if(e_l==Constant.CAP0 || e_l==Constant.CAP1) {
		//					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
		//					System.exit(1);
		//				}
		//				remains[index(e_l)]=false;
		//
		//				e_r=e2.regular[e_l];
		//				length+=2;
		//				if(e_r==Constant.CAP0) {System.out.println("0'0'"+length);break;}
		//				if(e_r==Constant.CAP1) {System.out.println("1'0'"+length);break;}
		//				remains[index(e_r)]=false;
		//			}
		//		}
		//
		//		for(Integer v:e2.cap1) {
		//			if(!remains[index(v)]) continue;
		//			e_r=v;
		//			remains[index(e_r)]=false;
		//
		//			int length=0;
		//			while(true) {
		//				e_l=e1.regular[e_r];
		//				if(e_l==Constant.CAP0 || e_l==Constant.CAP1) {
		//					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
		//					System.exit(1);
		//				}
		//				remains[index(e_l)]=false;
		//
		//				e_r=e2.regular[e_l];
		//				length+=2;
		//				if(e_r==Constant.CAP0) {
		//					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
		//					System.exit(1);
		//				}
		//				if(e_r==Constant.CAP1) {System.out.println("1'1'"+length);break;}
		//				remains[index(e_r)]=false;
		//			}
		//		}
		//
		//		for(int v:available) {
		//			if(!remains[index(v)]) continue;
		//			e_r=v;
		//			start=v;
		//			int length=0;
		//			while(true) {
		//				e_l=e2.regular[e_r];
		//				remains[index(e_l)]=false;
		//				remains[index(e_r)]=false;
		//				length+=2;
		//				e_r=e1.regular[e_l];
		//				if(e_r==start) {System.out.println("cycle"+length);break;}
		//			}
		//		}


		ArrayList<Integer> result=null;
		float best=-1;
		for(int left:available) {
			int right=e1.regular[left];
			if(right>=Constant.CAP1) continue;
			int l2=e2.regular[left], r2=e2.regular[right];
			if(l2>=Constant.CAP1 || r2>=Constant.CAP1) continue;
			int l3=e1.regular[l2], r3=e1.regular[r2];
			if(l3>=Constant.CAP1 || r3>=Constant.CAP1) continue;
			if(l3!=e2.regular[r3]) continue;
			float w=wa.regular[left][right]+wa.regular[l2][r2]+wa.regular[l3][r3];
			if(w>=0.5 && w>best) {
				best=w;
				result=new ArrayList<Integer>();
				result.add(left);result.add(right);result.add(l2);result.add(r2);result.add(l3);result.add(r3);
			}
		}

		for(int left:available) {
			int right=e2.regular[left];
			if(right>=Constant.CAP1) continue;
			int l2=e1.regular[left], r2=e1.regular[right];
			if(l2>=Constant.CAP1 || r2>=Constant.CAP1) continue;
			int l3=e2.regular[l2], r3=e2.regular[r2];
			if(l3>=Constant.CAP1 || r3>=Constant.CAP1) continue;
			if(l3!=e1.regular[r3]) continue;
			float w=wa.regular[left][right]+wa.regular[l2][r2]+wa.regular[l3][r3];
			if(w>=0.5 && w>best) {
				best=w;
				result=new ArrayList<Integer>();
				result.add(left);result.add(right);result.add(l2);result.add(r2);result.add(l3);result.add(r3);
			}
		}
		return result;
	}

	ArrayList<Integer> within_cycles (boolean opt) {
		// if optimal if true, the only the ones with non-negative weight are used
		// otherwise the ones with the greatest negative weights are used
		boolean optimal=opt;
		ArrayList<Integer> result=new ArrayList<Integer>();
		ArrayList<Integer> cycles=new ArrayList<Integer>();

		boolean remains[]=new boolean[available.size()];
		for(int i=0;i<available.size();i++) remains[i]=true;

		int e_l, e_r, start;
		for(Integer v:e1.cap0) {
			if(!remains[index(v)]) continue;
			cycles=new ArrayList<Integer>();
			e_r=v;
			remains[index(e_r)]=false;
			cycles.add(e_r);

			int length=0;
			while(true) {
				e_l=e2.regular[e_r];
				length++;
				if(e_l==Constant.CAP0) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path 00'",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
					break;
				}
				else if(e_l==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
					break;}
				remains[index(e_l)]=false;
				cycles.add(e_l);

				e_r=e1.regular[e_l];
				length++;
				if(e_r==Constant.CAP0) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
					break;}
				else if(e_r==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path 01",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path 01"));
					break;}
				remains[index(e_r)]=false;
				cycles.add(e_r);
			}
		}

		for(Integer v:e1.cap1) {
			if(!remains[index(v)]) continue;
			cycles=new ArrayList<Integer>();
			e_r=v;
			remains[index(e_r)]=false;
			cycles.add(e_r);

			int length=0;
			while(true) {
				e_l=e2.regular[e_r];
				length++;
				if(e_l==Constant.CAP0) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				else if(e_l==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				remains[index(e_l)]=false;
				cycles.add(e_l);

				e_r=e1.regular[e_l];
				length++;
				if(e_r==Constant.CAP0) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				else if(e_r==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				remains[index(e_r)]=false;
				cycles.add(e_r);
			}
		}

		for(Integer v:e2.cap0) {
			if(!remains[index(v)]) continue;
			cycles=new ArrayList<Integer>();
			e_r=v;
			remains[index(e_r)]=false;
			cycles.add(e_r);

			int length=0;
			while(true) {
				e_l=e1.regular[e_r];
				if(e_l==Constant.CAP0 || e_l==Constant.CAP1) {
					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
					System.exit(1);
				}
				remains[index(e_l)]=false;
				cycles.add(e_l);

				e_r=e2.regular[e_l];
				length+=2;
				if(e_r==Constant.CAP0) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				if(e_r==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path 0'1'",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path 0'1'"));
					break;}
				remains[index(e_r)]=false;
				cycles.add(e_r);
			}
		}

		for(Integer v:e2.cap1) {
			if(!remains[index(v)]) continue;
			cycles=new ArrayList<Integer>();
			e_r=v;
			remains[index(e_r)]=false;
			cycles.add(e_r);

			int length=0;
			while(true) {
				e_l=e1.regular[e_r];
				if(e_l==Constant.CAP0 || e_l==Constant.CAP1) {
					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
					System.exit(1);
				}
				remains[index(e_l)]=false;
				cycles.add(e_l);

				e_r=e2.regular[e_l];
				length+=2;
				if(e_r==Constant.CAP0) {
					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
					System.exit(1);
				}
				if(e_r==Constant.CAP1) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"path",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"path"));
					break;}
				remains[index(e_r)]=false;
				cycles.add(e_r);
			}
		}

		for(int v:available) {
			if(!remains[index(v)]) continue;
			cycles=new ArrayList<Integer>();
			e_r=v;
			start=v;
			int length=0;
			while(true) {
				e_l=e2.regular[e_r];
				remains[index(e_l)]=false;
				remains[index(e_r)]=false;
				cycles.add(e_r);
				cycles.add(e_l);
				length+=2;
				e_r=e1.regular[e_l];
				if(e_r==start) {
					ArrayList<Integer> choice=new ArrayList<Integer>();
					if(highest_weight(cycles,"cycle",choice)>=0) result.addAll(choice);
					else if(!optimal) result.addAll(choice);
//					result.addAll(highest_weight(cycles,"cycle"));
					break;}
			}
		}

		return result;
	}

	private double highest_weight(ArrayList<Integer> cycles, String type, ArrayList<Integer> result) {
		int cnt_cap=num_free_caps;
		boolean has_cap_in_final=false, has_cap;
		//		for(int v:cycles) 
		//			System.out.print(v+"\t");
		int left,right,lindex,rindex;
		ArrayList<Integer> best=new ArrayList<Integer>();
		double best_weight=-1; // weight per edge
		if(type.equals("cycle")) {
			int size=cycles.size();
			for(int i=0;i<size;i++) {
				//				System.out.println("considering edge "+i+" "+(i+1));
				lindex=i+1; rindex=i;
				float w=(float)-1;
				ArrayList<Integer> choice=new ArrayList<Integer>();
				for(int loop=1;loop<=size/2;loop++) {
					lindex--;rindex++;
					if(lindex<0) lindex+=size;
					if(rindex>=size) rindex-=size;
					left=cycles.get(lindex);
					right=cycles.get(rindex);
					choice.add(left);choice.add(right);
					if(e1.regular[left]==right || e2.regular[left]==right) {
						w+=wa.regular[left][right]+0.5;
					}
					else w+=(wa.regular[left][right]-0.5);
					double wpe= w;
					if(wpe>best_weight) {
						best=new ArrayList<Integer>(choice);
						best_weight=wpe;
					}
					//					System.out.println(loop+" "+w);
				}
				//				System.out.println();
			}
		}
		else if(type.equals("path 00'")) {
			cycles.add(Constant.CAP0);
			int size=cycles.size();
			out:for(int i=0;i<size;i++) {
				//				System.out.println("considering edge "+i+" "+(i+1));
				lindex=i+1; rindex=i;
				float w=(float)-1;
				ArrayList<Integer> choice=new ArrayList<Integer>();
				has_cap=false;
				for(int loop=1;loop<=size/2;loop++) {
					lindex--;rindex++;
					if(lindex<0) lindex+=size;
					if(rindex>=size) rindex-=size;
					left=cycles.get(lindex);
					right=cycles.get(rindex);
					if(left<Constant.CAP1 && right<Constant.CAP1) {
						choice.add(left);choice.add(right);
						if(e1.regular[left]==right || e2.regular[left]==right) {
							w+=wa.regular[left][right]+0.5;
						}
						else w+=(wa.regular[left][right]-0.5);
					}
					else if(left==Constant.CAP0 && cnt_cap>0) {
						choice.add(right);choice.add(left);
						has_cap=true;
						if(e1.cap0.contains((Integer) right) || e2.cap0.contains((Integer) right)) {
							w+=wa.cap0[right]+1-cap_thr2;
						}
						else w+=(wa.cap0[right]-0.5);
					}
					else if(right==Constant.CAP0 && cnt_cap>0) {
						choice.add(left);choice.add(right);
						has_cap=true;
						if(e1.cap0.contains((Integer) left) || e2.cap0.contains((Integer) left)) {
							w+=wa.cap0[left]+1-cap_thr2;
						}
						else w+=(wa.cap0[left]-0.5);
					}
					else continue out;

					double wpe= w;
					if(wpe>best_weight) {
						best=new ArrayList<Integer>(choice);
						best_weight=wpe;
						has_cap_in_final=has_cap;
					}
					//					System.out.println(loop+" "+w+" has_cap?"+has_cap_in_final);
				}
				//				System.out.println();
			}
			if(best_weight>=0 && has_cap_in_final) 
				cnt_cap--;
		}
		//					System.out.println(loop+" "+w);
		else if(type.equals("path 01") || type.equals("path 0'1'")) {
			int size=cycles.size();
			for(int i=0;i<size;i++) {
				//				System.out.println("considering edge "+i+" "+(i+1));
				lindex=i+1; rindex=i;
				float w=(float)-1;
				ArrayList<Integer> choice=new ArrayList<Integer>();
				for(int loop=1;loop<=size/2;loop++) {
					lindex--;rindex++;
					if(lindex<0) lindex+=size;
					if(rindex>=size) rindex-=size;
					left=cycles.get(lindex);
					right=cycles.get(rindex);
					choice.add(left);choice.add(right);
					if(e1.regular[left]==right || e2.regular[left]==right || lindex==0 && rindex==size-1 || lindex==size-1 && rindex==0) {
						w+=(wa.regular[left][right]+0.5);
					}
					else w+=(wa.regular[left][right]-0.5);


					double wpe= w;
					if(wpe>best_weight) {
						best=new ArrayList<Integer>(choice);
						best_weight=wpe;
					}
					//					System.out.println(loop+" "+w);
				}

				//				System.out.println();
			}
		}
		else {
			int size=cycles.size();
			for(int i=0;i<size;i++) {
				//				System.out.println("considering edge "+i+" "+(i+1));
				lindex=i+1; rindex=i;
				float w=(float)-1;
				ArrayList<Integer> choice=new ArrayList<Integer>();
				for(int loop=1;loop<=size/2;loop++) {
					lindex--;rindex++;
					if(lindex<0) continue;
					if(rindex>=size) continue;
					left=cycles.get(lindex);
					right=cycles.get(rindex);
					choice.add(left);choice.add(right);
					if(e1.regular[left]==right || e2.regular[left]==right) {
						w+=wa.regular[left][right]+0.5;
					}
					else w+=(wa.regular[left][right]-0.5);
					double wpe= w;
					if(wpe>best_weight) {
						best=new ArrayList<Integer>(choice);
						best_weight=wpe;
					}
					//					System.out.println(loop+" "+w);
				}
				//				System.out.println();
			}
		}

		//		System.out.println("in highest_weight"+best_weight);
		if(!best.isEmpty()) result.addAll(best);
		return best_weight;
	}

	void add_median_edges(ArrayList<Integer> result) {
//		for(int i=0;i<result.size();i++) {
//			for(int j=i+1;j<result.size();j++) 
//				if(result.get(i).equals(result.get(j)) && !result.get(i).equals((Integer) Constant.CAP0)) {
//					int i2,j2;
//					if(i%2==0) i2=i+1;
//					else i2=i-1;
//					if(j%2==0) j2=j+1;
//					else j2=j-1;
//					System.out.printf("old pair %d %d\tnew pair %d %d\n",result.get(i),result.get(i2),result.get(j),result.get(j2));
//				}
//		}
		for(Integer v:result) {
			if(v==Constant.CAP0) num_free_caps--;
			else available.remove((Integer) v); 
		}
		for(int i=0;i<result.size();i+=2) {
			int left=result.get(i), right=result.get(i+1);
			if(left>=Constant.NULL || right>=Constant.NULL ) 
				System.out.println();
			median.reg_adj[left]=right;
			if(right==Constant.CAP0) {
				median.cap_adj.add(left);
			}
			else median.reg_adj[right]=left;
		}
		size-=(result.size())/2;
		//		if(2*size!=available.size()+num_free_caps)
		//				System.out.println();
		e1.add_median_edges(result);
		e2.add_median_edges(result);
		wa.add_median_edges(result);

//		for(int i=0;i<2*num_gene;i++) 
//			if(!available.contains(i) && median.reg_adj[i]==Constant.NULL) 
//				System.out.println();
	}
}

class Edges {
	float weight;
	int num_gene;
	int regular[];
	ArrayList<Integer> cap0;
	ArrayList<Integer> cap1;

	public Edges (float w, Adjacency adj) {
		weight=w;
		num_gene=adj.num_gene;
		regular=new int[2*num_gene];
		cap0=new ArrayList<Integer>();
		cap1=new ArrayList<Integer>();


		for(int left=0;left<2*num_gene;left++) {
			int right=adj.reg_adj[left];
			regular[left]=right;
			if(right==Constant.CAP0) cap0.add(left);
			else if(right==Constant.CAP1) cap1.add(left);
			else if(right==Constant.NULL) ;
		}
	}

	public Edges (float w, AdjNode node) {
		this(w,node.adj);
	}

	public void print() {
		for(int left=0;left<2*num_gene;left++) {
			int right=regular[left];
			if(right<Constant.CAP1 && right>left) System.out.printf("%d %d\n",left, right);
		}
		for(int v:cap0) System.out.printf("%d cap0\n",v);
		for(int v:cap1) System.out.printf("%d cap1\n",v);
		System.out.println();
	}

	public void add_median_edges(ArrayList<Integer> zeros) {
		for(int i=0;i<zeros.size();i+=2) {
			int left=zeros.get(i), right=zeros.get(i+1);
			if(left<Constant.CAP1 && right<Constant.CAP1) {
				int l_c=regular[left];
				int r_c=regular[right];
				// a lot things to consider, depending whether l_c, r_c take values Constant.CAP0, Constant.CAP1, NULL

				if(l_c!=Constant.CAP0 && l_c!=Constant.CAP1 && l_c!=Constant.NULL) {
					if(r_c!=Constant.CAP0 && r_c!=Constant.CAP1 && r_c!=Constant.NULL) {
						regular[r_c]=l_c;
						regular[l_c]=r_c;
					}
					else if(r_c==Constant.CAP0) {
						cap0.remove((Integer)right);
						cap0.add(l_c);
						regular[l_c]=Constant.CAP0;
					}
					else if(r_c==Constant.CAP1) {
						cap1.remove((Integer)right);
						cap1.add(l_c);
						regular[l_c]=Constant.CAP1;
					}
					else if(r_c==Constant.NULL) {
						regular[l_c]=Constant.NULL;
					}
				}
				else if(l_c==Constant.CAP0) {
					cap0.remove((Integer)left);
					if(r_c!=Constant.CAP0 && r_c!=Constant.CAP1 && r_c!=Constant.NULL ) {
						cap0.add(r_c);
						regular[r_c]=Constant.CAP0;
					}
					else if(r_c==Constant.CAP0) {
						cap0.remove((Integer)right);
					}
					else if(r_c==Constant.CAP1) {
						cap1.remove((Integer)right);
					}
					else if(r_c==Constant.NULL) ;
				}
				else if(l_c==Constant.CAP1) {
					cap1.remove((Integer)left);
					if(r_c!=Constant.CAP0 && r_c!=Constant.CAP1 && r_c!=Constant.NULL ) {
						cap1.add(r_c);
						regular[r_c]=Constant.CAP1;
					}
					else if(r_c==Constant.CAP0) {
						cap0.remove((Integer)right);
					}
					else if(r_c==Constant.CAP1) {
						cap1.remove((Integer)right);
					}
					else if(r_c==Constant.NULL) ;
				}
				else if(l_c==Constant.NULL) {
					if(r_c!=Constant.CAP0 && r_c!=Constant.CAP1 && r_c!=Constant.NULL ) {
						regular[r_c]=Constant.NULL;
					}
					else if(r_c==Constant.CAP0) {
						cap0.remove((Integer)right);
					}
					else if(r_c==Constant.CAP1) {
						cap0.remove((Integer)right);
					}
					else if(r_c==Constant.NULL) ;
				}
				regular[left]=Constant.NULL;
				regular[right]=Constant.NULL;
			}
			else if(left<Constant.CAP1 && right==Constant.CAP0) {
				int l_c=regular[left];
				if(l_c!=Constant.CAP0 && l_c!=Constant.CAP1 && l_c!=Constant.NULL) {
					cap1.add(l_c);
					regular[l_c]=Constant.CAP1;
				}
				if(l_c==Constant.CAP0) {
					cap0.remove((Integer)left);
				}
				if(l_c==Constant.CAP1) {
					cap1.remove((Integer)left);
				}
				if(l_c==Constant.NULL) ;
				regular[left]=Constant.NULL;
			}
		}
	}
}

class WeightedAdjacency {
	int num_gene;
	float regular[][]=null;
	float cap0[];
	float cap1[];

	ArrayList<Edges> leaf_edges;

	public WeightedAdjacency() {
		leaf_edges=new ArrayList<Edges>();
	}

	public WeightedAdjacency(int num_g) {
		leaf_edges=new ArrayList<Edges>();
		num_gene=num_g;
		regular=new float[2*num_gene][2*num_gene];
		cap0=new float[2*num_gene];
		cap1=new float[2*num_gene];
		for(int i=0;i<2*num_gene;i++) {
			for(int j=0;j<2*num_gene;j++) 
				regular[i][j]=0;
			cap0[i]=0;
			cap1[i]=0;
		}
	}

	public void add_leaf(float w, AdjNode node) {
		if(regular==null) {
			num_gene=node.adj.num_gene;
			regular=new float[2*num_gene][2*num_gene];
			cap0=new float[2*num_gene];
			cap1=new float[2*num_gene];
			for(int i=0;i<2*num_gene;i++) {
				for(int j=0;j<2*num_gene;j++) 
					regular[i][j]=0;
				cap0[i]=0;
				cap1[i]=0;
			}
		}
		leaf_edges.add(new Edges(w,node));
	}

	public void update() {
		regular=new float[2*num_gene][2*num_gene];
		cap0=new float[2*num_gene];
		cap1=new float[2*num_gene];
		for(int i=0;i<2*num_gene;i++) {
			for(int j=0;j<2*num_gene;j++) 
				regular[i][j]=0;
			cap0[i]=0;
			cap1[i]=0;
		}

		for(Edges e:leaf_edges) {
			for(int left=0;left<2*num_gene;left++) {
				int right=e.regular[left];
				if(right==Constant.CAP0) cap0[left]+=e.weight;
				else if(right==Constant.CAP1) cap1[left]+=e.weight;
				else if(right==Constant.NULL || right<left) continue;
				else {
					regular[left][right]+=e.weight;
					regular[right][left]+=e.weight;
				}
			}
		}
	}

	public void add_median_edges(ArrayList<Integer> zeros) {
		for(int i=0;i<zeros.size();i+=2) {
			int left=zeros.get(i), right=zeros.get(i+1);
			for(Edges e: leaf_edges) {
				if(left<Constant.CAP1 && right<Constant.CAP1) {
					int l_c=e.regular[left];
					int r_c=e.regular[right];
					// a lot things to consider, depending whether l_c, r_c take values Constant.CAP0, Constant.CAP1, NULL

					if(l_c<Constant.CAP1 ) {
						regular[left][l_c]-=e.weight;
						regular[l_c][left]-=e.weight;
						if(r_c<Constant.CAP1 ) {
							e.regular[r_c]=l_c;
							e.regular[l_c]=r_c;
							regular[right][r_c]-=e.weight;
							regular[r_c][right]-=e.weight;
							regular[l_c][r_c]+=e.weight;
							regular[r_c][l_c]+=e.weight;
						}
						else if(r_c==Constant.CAP0) {
							e.cap0.remove((Integer)right);
							e.cap0.add(l_c);
							e.regular[l_c]=Constant.CAP0;
							cap0[right]-=e.weight;
							cap0[l_c]+=e.weight;
						}
						else if(r_c==Constant.CAP1) {
							e.cap1.remove((Integer)right);
							e.cap1.add(l_c);
							e.regular[l_c]=Constant.CAP1;
							cap1[right]-=e.weight;
							cap1[l_c]+=e.weight;
						}
						else if(r_c==Constant.NULL) {
							e.regular[l_c]=Constant.NULL;
						}
					}
					else if(l_c==Constant.CAP0) {
						e.cap0.remove((Integer)left);
						cap0[left]-=e.weight;
						if(r_c<Constant.CAP1 ) {
							e.cap0.add(r_c);
							e.regular[r_c]=Constant.CAP0;
							regular[right][r_c]-=e.weight;
							regular[r_c][right]-=e.weight;
							cap0[r_c]+=e.weight;
						}
						else if(r_c==Constant.CAP0) {
							e.cap0.remove((Integer)right);
							cap0[right]-=e.weight;
						}
						else if(r_c==Constant.CAP1) {
							e.cap1.remove((Integer)right);
							cap1[right]-=e.weight;
						}
						else if(r_c==Constant.NULL) ;
					}
					else if(l_c==Constant.CAP1) {
						e.cap1.remove((Integer)left);
						cap1[left]-=e.weight;
						if(r_c<Constant.CAP1 ) {
							e.cap1.add(r_c);
							e.regular[r_c]=Constant.CAP1;
							regular[right][r_c]-=e.weight;
							regular[r_c][right]-=e.weight;
							cap1[r_c]+=e.weight;
						}
						else if(r_c==Constant.CAP0) {
							e.cap0.remove((Integer)right);
							cap0[right]-=e.weight;
						}
						else if(r_c==Constant.CAP1) {
							e.cap1.remove((Integer)right);
							cap1[right]-=e.weight;
						}
						else if(r_c==Constant.NULL) ;
					}
					else if(l_c==Constant.NULL) {
						if(r_c!=Constant.CAP0 && r_c!=Constant.CAP1 && r_c!=Constant.NULL ) {
							e.regular[r_c]=Constant.NULL;
						}
						else if(r_c==Constant.CAP0) {
							e.cap0.remove((Integer)right);
						}
						else if(r_c==Constant.CAP1) {
							e.cap0.remove((Integer)right);
						}
						else if(r_c==Constant.NULL) ;
					}
					e.regular[left]=Constant.NULL;
					e.regular[right]=Constant.NULL;
				}
				else if(left<Constant.CAP1 && right==Constant.CAP0) {
					int l_c=e.regular[left];
					if(l_c<Constant.CAP1 ) {
						regular[left][l_c]-=e.weight;
						regular[l_c][left]-=e.weight;
						e.cap1.add(l_c);
						e.regular[l_c]=Constant.CAP1;
						cap0[l_c]+=e.weight;
					}
					if(l_c==Constant.CAP0) {
						e.cap0.remove((Integer)left);
						cap0[left]-=e.weight;
					}
					if(l_c==Constant.CAP1) {
						e.cap1.remove((Integer)left);
						cap1[left]-=e.weight;
					}
					if(l_c==Constant.NULL) ;
					e.regular[left]=Constant.NULL;
				}
			}
		}
	}

	public void print() {
		System.out.println();
		for(int i=0;i<2*num_gene;i++) {
			for(int j=i+1;j<2*num_gene;j++)
				if(regular[i][j]==1) System.out.printf("%4d %4d %8.4f\n", i,j,regular[i][j]);
		}
		for(int i=0;i<2*num_gene;i++) if(cap0[i]==1) System.out.printf("with caps %4d %8.4f\n", i,cap0[i]);
		for(int i=0;i<2*num_gene;i++) if(cap1[i]==1) System.out.printf("with cap0 %4d %8.4f\n", i,cap1[i]);
	}
}