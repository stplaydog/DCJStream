package gasts;

import java.util.Random;

class SimulatedGenome {
	final static short CAP=32767;
	int gene_number, chromosome_number;
	short genes[];
	
	public SimulatedGenome (int n, int xx) {
		gene_number=n;
		chromosome_number=xx;
		genes=new short[gene_number+chromosome_number];
		int gene_per_chromosome=gene_number/chromosome_number;
		int gene_in_last_chromosome=gene_number-(chromosome_number-1)*gene_per_chromosome;
		int count=0;
		for(int x=0;x<chromosome_number-1;x++) {
			for(int i=0;i<gene_per_chromosome;i++) 
				genes[count++]=(short)(x*gene_per_chromosome+i+1);
			genes[count++]=CAP;
		}
		for(int i=0;i<gene_in_last_chromosome;i++) 
			genes[count++]=(short)((chromosome_number-1)*gene_per_chromosome+i+1);
		genes[count]=CAP;
	}
	
	public SimulatedGenome (SimulatedGenome g) {
		gene_number=g.gene_number;
		chromosome_number=g.chromosome_number;
		genes=new short[gene_number+chromosome_number];
		for(int i=0;i<gene_number+chromosome_number;i++) 
			genes[i]=g.genes[i];
	}
	
	public void reverse(int r) {
		int count=0;
		while (count < r) {
			Random rndm = new Random();
			int left = rndm.nextInt(gene_number + chromosome_number-1);
			int right = rndm.nextInt(gene_number + chromosome_number-1);
			if (left > right) {
				int tmp = left;
				left = right;
				right = tmp;
			}
			if(genes[left]==CAP && genes[right]==CAP) continue;
			if((left==0 || genes[left-1]==CAP )&& genes[right+1]==CAP) continue;
			count++;
			short tmp[]=new short[gene_number+chromosome_number];
			for(int i=left;i<=right;i++) tmp[i]=genes[i];
			for(int i=left;i<=right;i++) {
				int opposite=right+left-i;
				if(tmp[opposite]==CAP) genes[i]=CAP;
				else genes[i]=(short)(tmp[opposite]*(-1));
			}
		}
	}
	
	public void tranpose(int r) {
		int count=0;
		out:while (count < r) {
			Random rndm = new Random();
			int left = rndm.nextInt(gene_number + chromosome_number);
			int right = rndm.nextInt(gene_number + chromosome_number);
			if(left==right) continue out;
			if (left > right) {
				int tmp = left;
				left = right;
				right = tmp;
			}
			if(left==0 && right==genes.length-1) continue;
			for(int i=left;i<=right;i++) 
				if(genes[i]==CAP) continue out;
			
			int third;
			int shift=right-left;
			if(rndm.nextFloat()<left/1.0/(genes.length-shift)) third=rndm.nextInt(left);
			else third=right+1+rndm.nextInt(genes.length-right-1);
			count++;
			
			
			short tmp[]=new short[shift];
			for(int i=left;i<right;i++) tmp[i-left]=genes[i];
			if(third>right) {
				for(int i=right;i<third;i++) genes[i-shift]=genes[i];
				for(int i=0;i<shift;i++) genes[third-shift+i]=tmp[i];
			}
			else {
				for(int i=left-1;i>=third;i--) genes[i+shift]=genes[i];
				for(int i=0;i<shift;i++) genes[third+i]=tmp[i];
			}
//			System.out.println("pause");
		}
	}
	
	public String get_string(String name) {
		StringBuilder sb =new StringBuilder();
		sb.append(">");sb.append(name);sb.append("\n");
		for(int i=0;i<gene_number + chromosome_number;i++) {
			if(genes[i]==CAP) sb.append("\n"); // CHANGE usually there is a $ at the end
			else {sb.append(genes[i]);sb.append(" ");}
		}
		return sb.toString();
	}
}
