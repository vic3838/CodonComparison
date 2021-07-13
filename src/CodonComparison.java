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
		//System.out.println(path);
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
    