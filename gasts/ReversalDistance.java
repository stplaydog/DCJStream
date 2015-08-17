package gasts;


import java.util.Random;
import java.util.Scanner;
import java.io.*;

public class ReversalDistance {
	static public int reversalDistance (SimuNode n1, SimuNode n2) throws Exception {
		Random rd=new Random();
		String file="tmp"+rd.nextInt();
		PrintWriter pw=new PrintWriter(file);
		pw.append(n1.genome.get_string(n1.name));
		pw.append(n2.genome.get_string(n2.name));
		pw.append(n2.genome.get_string("A4"));
		pw.close();
		
		Process process = Runtime.getRuntime().exec("./grimm -C -f "+file);
		BufferedReader input =new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line;
		Scanner sc;
		int d=0;
		while ((line = input.readLine()) != null) {
			if(line.contains("A4")) {
				sc=new Scanner(line);
				sc.skip("A4");
				d= sc.nextInt();
				break;
			}
		}
		input.close();
		File f=new File(file);
		f.delete();
		return d;
	}
	
	static public int reversalDistance (AdjNode n1, AdjNode n2, String option) throws Exception {
		Random rd=new Random();
		String file="tmp"+rd.nextInt();
		PrintWriter pw=new PrintWriter(file);
		pw.append(n1.adj.toGenome().toString());
		pw.append(n2.adj.toGenome().toString());
		Adjacency n3=new Adjacency(n2.adj);
		n3.name="A4";
		pw.append(n3.toGenome().toString());
		pw.close();
		
		Process process = Runtime.getRuntime().exec("./grimm "+option +" -f "+file);
		BufferedReader input =new BufferedReader(new InputStreamReader(process.getInputStream()));
		String line;
		Scanner sc;
		int d=0;
		while ((line = input.readLine()) != null) {
			if(line.contains("A4")) {
				sc=new Scanner(line);
				sc.skip("A4");
				d= sc.nextInt();
				break;
			}
		}
		input.close();
		File f=new File(file);
		f.delete();
		return d;
	}
}
