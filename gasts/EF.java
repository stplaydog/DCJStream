package gasts;

import java.util.ArrayList;
public class EF {
	ArrayList<ArrayList<Integer>> black;
	String info;
	public EF(){
		black=new ArrayList<ArrayList<Integer>>();
	}
	public EF(EF ef){
		black=new ArrayList<ArrayList<Integer>>();
		for(ArrayList<Integer> m:ef.black) black.add(new ArrayList<Integer>(m));
		info=ef.info;
	}
}
