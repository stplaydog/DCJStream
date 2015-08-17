package gasts;

import java.util.*;

public class Solution {
	int cycle_number;
	int heuristic;
	int upper_bound;
	int lower_bound;
	//long size;
	long time;
	long h_time;
	boolean finished;
	ArrayList<Integer> median;
	
	public Solution() {
		finished=true;
		median=new ArrayList<Integer>();
	}
	
	public String toString() {
		String s="C:"+cycle_number+"\tU:"+upper_bound+"\tL:"+lower_bound+"\tH:"+heuristic+"\tT:"+time+"\tHT:"+h_time;
		if(finished) 
			return "Solved.\t"+s;
		else
			return "Unsolved.\t"+s;
	}
}
