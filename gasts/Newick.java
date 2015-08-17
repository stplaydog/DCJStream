package gasts;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;

public class Newick<T extends Node> {
	final int Parent=Constant.Parent;
	int NODE_ID=-1;
	String format=null;
	T root=null;

	
	
	public Newick(String input) {
		format=input;
	}
	
	public Newick(File file) throws Exception {
		BufferedReader br=new BufferedReader(new FileReader(file));
		format=br.readLine();
	}

	public  T readTree(Class<T> _class) throws Exception { // 
		// the input string does not contain a ";" at the end
		T current_node=null, next_node, previous=null;
		int name_begins=-1;
		for(int index=0;index<format.length();index++) {
			if(format.charAt(index)=='(')  {
				name_begins=index+1;
				NODE_ID++;
				next_node=(T) _class.newInstance();
				next_node.set_name("_"+NODE_ID);
				if(root==null) {
					root=next_node;
					current_node=next_node;
				}
				else {
					current_node.add_child(next_node);
					current_node=next_node;
				}
			}
			else if(format.charAt(index)==',') {
				if(index==name_begins) ;
				else if(format.charAt(name_begins)==':') { // edge_distance
					previous.edge_length=Double.parseDouble(format.substring(name_begins+1,index));
				}
				else  {
					next_node=(T) _class.newInstance();
					String code=format.substring(name_begins,index);
					int colon=code.indexOf(':');
					if(colon==-1) {
						next_node.set_name(code);
						next_node.edge_length=1;
					}
					else {
						next_node.set_name(format.substring(name_begins,colon+name_begins));
						next_node.edge_length=Double.parseDouble(format.substring(colon+1+name_begins,index));
					}
					next_node.setLeaf();
					current_node.add_child(next_node);
				}
				name_begins=index+1;
			}
			else if(format.charAt(index)==')') {
				if(index==name_begins) ; // nothing happens
				else if(format.charAt(name_begins)==':') { // edge_distance
					previous.edge_length=Double.parseDouble(format.substring(name_begins+1,index));
				}
				else { // adding a leaf node
					next_node=(T) _class.newInstance();
					String code=format.substring(name_begins,index);
					int colon=code.indexOf(':');
					if(colon==-1) {
						next_node.set_name(code);
						next_node.edge_length=1;
					}
					else {
						next_node.set_name(format.substring(name_begins,colon+name_begins));
						next_node.edge_length=Double.parseDouble(format.substring(colon+1+name_begins,index));
					}
					next_node.setLeaf();
					current_node.add_child(next_node);
				}
				previous=current_node;
				current_node=(T) current_node.links.get(Parent);
				name_begins=index+1;
			}
			else if(format.charAt(index)==';') break;
		}

		return root;
	}
		
	public String get_Newick(Node node, boolean has_length) {
		String f="";
		if(node.links.size()==1) {
			if(has_length) f=f+node.name+":"+node.edge_length;
			else f+=node.name;
		}
		else {
			f+="(";
			for(int i=1;i<node.links.size();i++) {
				f+=get_Newick(node.links.get(i),has_length);
				if(i<node.links.size()-1) f+=",";
			}
			if(has_length) f=f+"):"+node.edge_length;
			else f+=")";
		}
		return f;
	}
	
	public String get_Newick(boolean has_length) {
		if(root.isLeaf) {
			String f; 
			if(has_length) f="("+root.name+":"+root.edge_length+",";
			else f="("+root.name+",";
			for(int i=1;i<root.links.size();i++) {
				f+=get_Newick(root.links.get(i),has_length);
				if(i<root.links.size()-1) f+=",";
			}
			f+=")";
			return f;
		}
		else return get_Newick(root, has_length);
	}
	
	
//	public T find_node_by_name(String name) {
//		return find_node_by_name(root, name);
//	}
//	
//	public T find_node_by_name(T node, String name) {
//		System.out.println(node.name+" "+name);
//		if(node.name.equals(name)) return node;
//		
//		if(node.num_ch==0) return null;
//		else if(node.num_ch==1) return find_node_by_name(node.link[0],name);
//		else {
//			Node result=find_node_by_name(node.link[0],name);
//			if(result!=null) return result;
//			else return find_node_by_name(node.link[1],name);
//		}
//	}
	
	public T reroot () {
		T left=(T) root.links.get(1);
		T right=(T) root.links.get(2);
		right.cut();left.cut();
		right.edge_length+=left.edge_length;
		left.edge_length=0;
		
		while(left!=null) {
			T new_left;
			if(left.isLeaf) new_left=null;
			else new_left=(T) left.links.get(1);
			
			if(new_left!=null) {new_left.cut();left.edge_length=new_left.edge_length;new_left.edge_length=0;}
			left.add_child(right);
			right=left;
			left=new_left;
		}
		
		root=right;
		return root;
	}
	
	public void unroot() {
		T left=(T) root.links.get(1);
		T right=(T) root.links.get(2);
		right.cut();left.cut();
		right.edge_length+=left.edge_length;
		left.edge_length=0;
		
		while(left!=null) {
			T new_left;
			if(left.isLeaf) new_left=null;
			else new_left=(T) left.links.get(1);
			
			if(new_left!=null) {new_left.cut();left.edge_length=new_left.edge_length;new_left.edge_length=0;}
			left.add_child(right);
			right=left;
			left=new_left;
		}
		
		root=right;
	}
	
	public static void main(String[] args) throws Exception {
		Newick<Node> nw=new Newick<Node>(new File("tree"));
		nw.readTree(Node.class);
		System.out.println(nw.get_Newick(true));
		nw.unroot();
		System.out.println(nw.get_Newick(true));
	}
}



