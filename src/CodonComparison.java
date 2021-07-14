import java.util.*;
import java.io.*;

//@Author Vijay Kiran Cherupally
//god please help me
public class CodonComparison {
	
	static final String START = "ATG";
	static final String STOP1 = "TGA";
	static final String STOP2 = "TAA";
	static final String STOP3 = "TAG";
	

	public static void main(String[] args) throws IOException {
		//this method will print out the LCA sequence reconstruction from every tool
		getEmergingORFs();
		
		
	}
	
	
	public static void getEmergingORFs() throws IOException
	{
		String path = "D:\\College\\TECBio\\CodonComparison";	//create a string with the path of the current directory, will be used and edited
		File dir = new File(path);								//creating the file of the current directory
		System.out.println(dir.getAbsolutePath());				//testing shit
		
		File [] emergingORFs = dir.listFiles(new MyFileNameFilter()); 	//get every item in this directory except for the git, bin, and src directories
		//going through the emergingORFs, YAL, YBr, ...
		for(File x : emergingORFs)
		{		
			if(x.isDirectory())
			{
				//return an arraylist of the 8 reconstructed sequences and print it
				getReconstructions(x.getAbsolutePath());
				//System.out.println(x.toString());
			}
			
		}
	}
	
	//this method gets and prints the 8 reconstructions of each locus, might make it a return method later
	public static void getReconstructions(String path) throws IOException
	{
		String newPath = "";	//will be used to pass in the path to the subdirectories
		//in each locus, there are the eight reconstructions, write 8 helper methods to return the LCA string	
		System.out.println(path);
		File dir = new File(path);
		String [] reconDirs = dir.list(new MyDirectoryFilter());			//filters out any files that have extensions, leaving only directories
		String [] recons = new String[8];
		for(int i = 0; i < 8; i++)
		{
			switch(i) {
			case 0: recons[i] = fastML_free_marg(path + "\\" + reconDirs[0]); 
					System.out.println("fastML_free_marg: " + recons[i]);
					break;
			case 1: recons[i] = fastML_free_joint(path + "\\" + reconDirs[0]); 
					System.out.println("fastML_free_joint: " + recons[i]);
					break;
			case 2: recons[i] = fastML_sp_marg(path + "\\" + reconDirs[1]); 
					System.out.println("fastML_sp_marg: " + recons[i]);
					break;
			case 3: recons[i] = fastML_sp_joint(path + "\\" + reconDirs[1]); 
					System.out.println("fastML_sp_joint: " + recons[i]);
					break;
			case 4: recons[i] = prank_free(path + "\\" + reconDirs[2]); 
					System.out.println("prank_free: " + recons[i]);
					break;
			case 5: recons[i] = prank_sp(path + "\\" + reconDirs[3]); 
					System.out.println("prank_joint: " + recons[i]); 
					break;
			case 6: recons[i] = prequel_free(path + "\\" + reconDirs[4]); 
					System.out.println("prequel_free: " + recons[i]);
					break;
			case 7: recons[i] = prequel_sp(path + "\\" + reconDirs[5]); 
					System.out.println("prequel_sp: " + recons[i]); 
					break;
			default: System.out.println("ERROR");
			}
		}
		
		/*
		for(String x : recons)
		{
			System.out.println(x);
		}
		*/
		
		
	}

	//helper methods
	private static String prequel_free(String path)
	{
		return path;
	}
	private static String prequel_sp(String path)	throws IOException
	{
		path += "\\ORF_alignement.Seub-Sarb.fa";		//create the path for the file we need to open	
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		while(fileScan.hasNextLine()){
			fileLine = fileScan.nextLine();		
			if(fileLine.charAt(0) != '>')
			{
				seq += fileLine;
			}
		}
		
		return seq;
	}
	private static String fastML_free_marg(String path) throws IOException
	{
		path += "\\seq.marginal.txt";
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		boolean foundSeq = false;
		while(fileScan.hasNextLine()) {
			fileLine = fileScan.nextLine();
			if(foundSeq) {					//if foundSeq = true, that means we have found the sequence and this next line is the one we need
				seq = fileLine;
				break;
			}
			else if(fileLine.contains(">N1")) {
				foundSeq = true;
			}			
		}		
		return seq;
	}
	private static String fastML_free_joint(String path) throws IOException
	{
		path += "\\seq.joint.txt";
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		boolean foundSeq = false;
		while(fileScan.hasNextLine()) {
			fileLine = fileScan.nextLine();
			if(foundSeq) {					//if foundSeq = true, that means we have found the sequence and this next line is the one we need
				seq = fileLine;
				break;
			}
			else if(fileLine.contains(">N1")) {
				foundSeq = true;
			}			
		}		
		return seq;
	}
	private static String fastML_sp_marg(String path) throws IOException
	{ 
		path += "\\seq.marginal.txt";
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		boolean foundSeq = false;
		while(fileScan.hasNextLine()) {
			fileLine = fileScan.nextLine();
			if(foundSeq) {					//if foundSeq = true, that means we have found the sequence and this next line is the one we need
				seq = fileLine;
				break;
			}
			else if(fileLine.contains(">N1")) {
				foundSeq = true;
			}			
		}		
		return seq;
	}
	private static String fastML_sp_joint(String path) throws IOException
	{
		path += "\\seq.joint.txt";
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		boolean foundSeq = false;
		while(fileScan.hasNextLine()) {
			fileLine = fileScan.nextLine();
			if(foundSeq) {					//if foundSeq = true, that means we have found the sequence and this next line is the one we need
				seq = fileLine;
				break;
			}
			else if(fileLine.contains(">N1")) {
				foundSeq = true;
			}			
		}		
		return seq;
	}
	
	//need to get the largest identified node
	private static String prank_free(String path) throws IOException{
		path += "\\ORF_alignment.anc.fas";
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		//arrays to hold the sequences
		String[][] seqsAndNodes = new String[2][10];
		
		boolean storeSeq = false;
		int col = 0;
		//need to read through the file, find the highest numbered node, and then choose the sequence related to that...
		//probably add all the reconstructed seqs to an array and then choose the seq after so i dont have to read the file again
		while(fileScan.hasNextLine()){
			fileLine = fileScan.nextLine();
			if(storeSeq) {
				seqsAndNodes[1][col] = fileLine;
				col++;
				storeSeq = false;
			}
			
			else if(fileLine.contains(">#")) {
				seqsAndNodes[0][col] = fileLine;
				storeSeq = true;
			}				
		}
		
		
		
		return path;
		
	}
	private static String prank_sp(String path)
	{
		return path;
	}
	
	
	
	
	
	//filename filter override method to make sure we ignore certain directories, specifically .git, bin, and arc
	//FileNameFilter implementation
	public static class MyFileNameFilter implements FilenameFilter {

		@Override
		public boolean accept(File dir, String name) {
			
			if(name.equals(".git") || name.equals("bin") || name.equals("src"))
			{
				return false;
			}
			
			else {
				return true;
			}
			
		}
	}
	
	public static class MyDirectoryFilter implements FilenameFilter {

		@Override
		public boolean accept(File dir, String name) {
			
			if(name.contains("."))
			{
				return false;
			}
			
			else {
				return true;
			}
			
		}

	}

	//method to check if a codon is a stop codon
	public static boolean isStop(String codon)
	{
		if(codon.equals(STOP1) || codon.equals(STOP2) || codon.equals(STOP3))
		{
			return true;
		}
		else 
		{
			return false;
		}
	
	}
	
}
    