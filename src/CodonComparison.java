import java.util.*;
import java.io.*;

//@Author Vijay Kiran Cherupally
//god please help me
public class CodonComparison {
	
	static final String START = "ATG";
	static final String STOP1 = "TGA";
	static final String STOP2 = "TAA";
	static final String STOP3 = "TAG";
	

	public static void main(String[] args) {
		//this method will print out the LCA sequence reconstruction from every tool
		getEmergingORFs();
		
		
	}
	
	
	public static void getEmergingORFs()
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
	public static void getReconstructions(String path)
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
		File file = new File(path);						//file we need to enter 
		FileWriter fwrite = new FileWriter(file);		
		PrintWriter writer = new PrintWriter(fwrite);
		
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		boolean newSequence = false;
		while(fileScan.hasNextLine()){
			fileLine = fileScan.nextLine();				//this is the current line
			if(fileLine.charAt(0) == '>')
			{
				newSequence = true;
			}
			//if it is a sequence and we have just started a new sequence, add the line to fastasequence
			if(fileLine.charAt(0) != '>' && newSequence)
			{
				fileLine = fileLine.toUpperCase();
				newSequence = false;
			}
		}
		
		
		writer.close();
		
		return path;
	}
	private static String fastML_free_marg(String path)
	{
		return path;
	}
	private static String fastML_free_joint(String path)
	{
		return path;
	}
	private static String fastML_sp_marg(String path)
	{
		return path;
	}
	private static String fastML_sp_joint(String path)
	{
		return path;
	}
	private static String prank_free(String path)
	{
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
    