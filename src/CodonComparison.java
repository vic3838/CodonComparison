import java.util.*;
import java.io.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.template.ProfileProfileAligner; 
import org.biojava.nbio.alignment.Alignments.ProfileProfileAlignerType;
import org.biojava.nbio.core.alignment.SimpleProfile;
import org.biojava.nbio.core.alignment.template.Profile;
import org.biojava.nbio.core.alignment.template.AlignedSequence;
import org.biojava.nbio.core.sequence.DNASequence; 
import org.biojava.nbio.core.sequence.MultipleSequenceAlignment;
import org.biojava.nbio.core.sequence.compound.NucleotideCompound; 
import org.biojava.nbio.core.sequence.io.FastaReaderHelper; 
import org.biojava.nbio.core.util.ConcurrencyTools;

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
		align();
		
		
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
				System.out.println();
				//System.out.println(x.toString());
			}
			
		}
	}
	
	//this method gets and prints the 8 reconstructions of each locus, might make it a return method later
	public static void getReconstructions(String path) throws IOException
	{
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
					System.out.println("prank_sp: " + recons[i]); 
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
		
		//this is to get the name of the highest node, which is the LCA, and find its position in the 2d array
		//once finding this index, we can select the corresponding string
		String loc = ">#" + col + "#";
		int index = -1;
		for(int i = 0; i < 10; i++)
		{
			if(seqsAndNodes[0][i].equals(loc)) {
				index = i;
				break;
			}
		}
		
		seq = seqsAndNodes[1][index];
		return seq;		
	}
	
	
	private static String prank_sp(String path) throws IOException
	{
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
		
		//this is to get the name of the highest node, which is the LCA, and find its position in the 2d array
		//once finding this index, we can select the corresponding string
		String loc = ">#" + col + "#";
		int index = -1;
		for(int i = 0; i < 10; i++)
		{
			if(seqsAndNodes[0][i].equals(loc)) {
				index = i;
				break;
			}
		}
		
		seq = seqsAndNodes[1][index];
		return seq;	
	}
	
	//this method will lead the alignments of the whatever sequences we give it
	public static void align() 
	{
		try {
		//ybr seqs for test
		//fastmlfreemarg
		DNASequence sequence1 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCCATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCAATTGGTTTTGCAGACGGTATTTTTCAATTCCTTTTTTAGGTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCGTATCTTAA");
		//fastML_free_joint
		DNASequence sequence2 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCCATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTTGGTTTTGCAGACGGTATTTTTCAATTCCTTTTTTGGGTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCGTATCTTAA");
		//fastmlspmarg
		DNASequence sequence3 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCTATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTCGGTTTTGCAGACGGTATTTTTTAATTCCTTTTTTAGGTTTTGTTTCGTCTTTCCTTTTTTTTATTGTTTTCGTATCTTAA");
		//fastML_sp_joint
		DNASequence sequence4 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCTATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTCGGTTTTGCAGACGGTATTTTTTAATTCCTTTTTTGGGTTTTGTTTCGTCTTTCCTTTTTTTTATTGTTTTCGTATCTTAA");
		//prank_free
		DNASequence sequence5 = new DNASequence("ATGTCCCATCTCTATATATATTCATTGACGG--CATTCTATTTTCTTGGT----ATTGAAATGAGCGTTTTTTATTATTACAATTGGTTTTGCAGACGGCACTTTCC--TTTCCTTTT---GTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCATATCTTAA");	
		//pranksp
		DNASequence sequence6 = new DNASequence("GAGGCTCCTGTCTATATATACCATTTGACATTCTATTCTATTCTACAACT----GTTAAAGTGAGTGGTTTTTATTGTAGAAGCTGGTTTTAAGGGTGGTATTTTTT--CTCCTTCCCTGGGTTTTGT------------ATATTTTATTGTTCTCCTATCTTAA");
		//prequel_sp
		DNASequence sequence7 = new DNASequence("GTGTTCCCTGTCTATATATATTCATTGACGGTATTTGTTTCTTTTGGTGTTGAAGTGAGCGTTTTTTATTATTGCACTTGGTTTTCCAGACGGTATTTTAATTTTATTTTTAAGTCTTGTTCTTCTTGCCTTTTTTTCTTTGTTCTGGTATCTTAA");	
		
		
		List<DNASequence> lst = new ArrayList<DNASequence>();
		lst.add(sequence1);
		lst.add(sequence2);
		lst.add(sequence3);
		lst.add(sequence4);
		lst.add(sequence5);
		lst.add(sequence6);
		lst.add(sequence7);
		
		Profile<DNASequence, NucleotideCompound> results = Alignments.getMultipleSequenceAlignment(lst);
		
		List<AlignedSequence<DNASequence, NucleotideCompound>> huuuuh = results.getAlignedSequences();
		for(AlignedSequence<DNASequence, NucleotideCompound> x : huuuuh)
		{
			System.out.println(x);
		}
		
		}
		catch(Exception e){
			e.printStackTrace();
		}
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
    