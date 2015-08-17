package gasts;


import java.util.Vector;

public class Node  {
	public double edge_length;
	public double start_time;
	public static int Parent=Constant.Parent;
	public String name;
	public boolean isLeaf=false;
	public Vector<Node> links =new Vector<Node>();
	
	public Node() {
		links.add(Parent,null); // prevent assigning any value here.
		edge_length=1;
	}
		
	public void setLeaf() {
		isLeaf=true;
	}
	
	public void set_name(String nm) {
		name=nm;
	}
	
	public void add_child(Node node) {
		node.links.set(Parent, this);
		links.add(node);
	}
	
	public void cut() {
		links.get(Parent).links.remove(this);
		links.set(Parent, null);
	}
	
	public Node unroot() {
		if(links.get(Parent)!=null) return null;
		Node left=links.get(1);
		Node right=links.get(2);
		right.cut();left.cut();
		right.edge_length+=left.edge_length;
		left.edge_length=0;
		
		while(left!=null) {
			Node new_left;
			if(left.isLeaf) new_left=null;
			else new_left=left.links.get(1);
			
			if(new_left!=null) {new_left.cut();left.edge_length=new_left.edge_length;new_left.edge_length=0;}
			left.add_child(right);
			right=left;
			left=new_left;
		}
		
		return right;
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
		if(links.get(Parent)!=null) return null;
		if(isLeaf) {
			String f; 
			if(has_length) f="("+name+":"+edge_length+",";
			else f="("+name+",";
			for(int i=1;i<links.size();i++) {
				f+=get_Newick(links.get(i),has_length);
				if(i<links.size()-1) f+=",";
			}
			f+=");";
			return f;
		}
		else return get_Newick(this, has_length);
	}
}
