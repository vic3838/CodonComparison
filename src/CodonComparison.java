import java.util.*;
import java.io.*;
import org.biojava.nbio.alignment.Alignments;
import org.biojava.nbio.alignment.Alignments.PairwiseSequenceAlignerType;
import org.biojava.nbio.core.alignment.matrices.SubstitutionMatrixHelper;
import org.biojava.nbio.core.alignment.template.SubstitutionMatrix;
import org.biojava.nbio.alignment.SimpleGapPenalty;
import org.biojava.nbio.core.alignment.template.SequencePair;
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
//		List<AlignedSequence<DNASequence, NucleotideCompound>> alignedSeqs = align();
//		for(AlignedSequence<DNASequence, NucleotideCompound> x : alignedSeqs)
//		{
//			System.out.println(x.toString());
//		}
		
		
		//TODO: make the void methods getEmergingORFs and align return actual values that can be operated on in the main method, will need to loop
		//		might have to call align() from getEmergingORFs()
		//now that i have the aligned sequences, the next basic step is to get the # of aligned codons and output them to a table
		
		
		
		
	}
	
	//take 2 sequences, reconstruction and extant sequence, and find the number of conserved codons 
	public static int[] numberOfAlignedCodons(String recon, String extant)
	{
		//since the strings are aligned, we can read codon by codon. worry about finding an ORF later
		StringBuilder rec = new StringBuilder(recon);
		StringBuilder ext = new StringBuilder(extant);
		int codons = 0;
		int alignedCodons = 0;
		int length = recon.length();
		String codon1 = null;
		String codon2 = null;
		for(int i = 0; i < length; i += 3) {
			
			if(i + 3 <= length) {
				codon1 = rec.substring(i, i + 3);
				codon2 = ext.substring(i, i + 3);
				
				if(codon1.equals(codon2)) {
					alignedCodons++;
				}
				
				codons++;
			}
		}
		if(codons != length/3) System.out.println("ERROR");	
		return new int[] {alignedCodons, codons};
	}
	
	
	public static void getEmergingORFs() throws IOException
	{
		String path = "D:\\College\\TECBio\\CodonComparison";	//create a string with the path of the current directory, will be used and edited
		File dir = new File(path);								//creating the file of the current directory
		System.out.println(dir.getAbsolutePath());				//testing shit
		
		File [] emergingORFs = dir.listFiles(new MyFileNameFilter()); 	//get every item in this directory except for the git, bin, and src directories
		String[] recons;
		String recon;													//holder for the recon
		String extant;													//holder for the extant seq
		AlignedSequence<DNASequence, NucleotideCompound> rec;			//holders for the aligned seqs
		AlignedSequence<DNASequence, NucleotideCompound> ext;
		int[] results;
		int count;
		String locus;													 //string that holds the ORF name (YBR, YAL, etc))
		//going through the emergingORFs, YAL, YBr, ...
		for(File x : emergingORFs)
		{		
			String str = x.getAbsolutePath();				//path to the ORF folder
			if(x.isDirectory())
			{				
				count = 0;				//count will be the number passed into outTable(), corresponds to the method being printed
				//return an array of the 8 reconstructed sequences and print it
				recons = getReconstructions(str);
				locus = str.substring(str.lastIndexOf("\\"));
				extant = str + locus + "_ali.fa";	//set extant to the path to the file containing the extant Scer sequence
				extant = getExtantSequence(extant);									//redefine extant to the extant sequence itself
				System.out.println("Extant seq: " + extant);
				
				//output table file
				File newFile = new File(str + "\\" + locus + "_alignedCodonsByLongestCodingSequence.txt");
				FileWriter fwrite = new FileWriter(newFile);
				PrintWriter writer = new PrintWriter(fwrite);
				//set up the file
				writer.println(locus + " aligned codons");
				writer.println("____________________");
				System.out.println();
				
				//pairwise alignment of each reconstruction with the extant sequence
				for(String y: recons) {
					//get the longest ORF/coding sequence from the reconstruction
					//find all the start indices in the reconstruction
					ArrayList<Integer> startIndices = findAllStart(y);
					//run findCodingSequence on all these start codons
					ArrayList<String> codingSequences = new ArrayList<String>();
					//run through the indexes of the start codons, find the ORF corresponding to it
					for(int j = 0; j < startIndices.size(); j++)
					{
						recon = findORF(y, startIndices.get(j));		//here recon is used as a holder for the ORF being found
						codingSequences.add(recon);
					}					
					//now, codingSequences should be full of the potential coding sequences
					//just select the longest of them and use that as our coding sequence
					Collections.sort(codingSequences, Comparator.comparing(String :: length));
					//use recon again to get the largest coding sequence from the reconstruction and use that in our alignment
					
					
					//if the table has an entry, we probably found a coding seq, but it could be "BAD"
					if(codingSequences.size() != 0) {
						recon = codingSequences.get(codingSequences.size() - 1);	
						//"BAD" means no coding sequence was found
						if(recon.equals("BAD")) {
							recon = "NA";
							results = new int[] {-1,-1};
							outTable(writer, results, count);
							count++;
						}
						//otherwise, we can proceed with the alignment
						else { 
							SequencePair<DNASequence, NucleotideCompound> alignment = align(recon, extant);	//align the sequences
							rec = alignment.getQuery();
							ext = alignment.getTarget();
							recon = rec.toString();						//string representation of the reconstruction alignment
							extant = ext.toString();					//string representation of the extant alignment
							results = numberOfAlignedCodons(recon, extant);
							outTable(writer, results, count);
							count++;
						}

					}
					//if the table is empty, no ORF was found
					else {
						recon = "NA";
						results = new int[] {-1,-1};
						outTable(writer, results, count);
						count++;
					}


		
					
					System.out.println(recon);
					System.out.println(extant);
					System.out.println(results[0] + "/" + results[1]);
					System.out.println();
					
				}
				
				
//				List<AlignedSequence<DNASequence, NucleotideCompound>> alignedSeqs = align(recons, extant);
//				for(AlignedSequence<DNASequence, NucleotideCompound> y : alignedSeqs)
//				{
//					System.out.println(y.toString());
//					
//				}
//
//				//use list iterator to save on runtime, get the last item of the list (extant sequence first)
//				extant = alignedSeqs.get(alignedSeqs.size() - 1).toString();
//				System.out.println("TEST: " + extant);
//				Iterator<AlignedSequence<DNASequence, NucleotideCompound>> it = alignedSeqs.iterator();
//				//while
//				
//				
//				System.out.println();

				writer.close();
			}
			
		}
	}
	
	//this method gets and prints the 8 reconstructions of each locus, might make it a return method later
	public static String[] getReconstructions(String path) throws IOException
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
			//TESTING PURPOSES ONLY: REVERT TO NORMAL AFTER TESTING, SHOULD BE PREQUEL FREE
			case 6: recons[i] = prequel_sp(path + "\\" + reconDirs[5]); 
					System.out.println("prequel_sp: " + recons[i]); 
					
					//recons[i] = prequel_free(path + "\\" + reconDirs[4]); 
					//System.out.println("prequel_free: " + recons[i]);
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
		
		return recons;
		
	}
	
	//
	public static String getExtantSequence(String path) throws IOException
	{
		Scanner fileScan = new Scanner(new FileInputStream(path));
		String fileLine = "";
		String seq = "";
		
		boolean foundSeq = false;
		while(fileScan.hasNextLine()){
			fileLine = fileScan.nextLine();		
			if(fileLine.charAt(0) == '>' && foundSeq) {
				break;
			}
			if(fileLine.equals(">Scer")){
				foundSeq = true;
			}
			if(fileLine.charAt(0) != '>' && foundSeq) {
				fileLine = fileLine.toUpperCase();
				seq += fileLine;
			}			
		}
		return seq;		
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
	public static SequencePair<DNASequence, NucleotideCompound> align(String recons, String extant) 
	{
		SequencePair<DNASequence, NucleotideCompound> psa = null;
		try {
			DNASequence rec = new DNASequence(recons);
			DNASequence ext = new DNASequence(extant);
				
			SubstitutionMatrix<NucleotideCompound> matrix = SubstitutionMatrixHelper.getNuc4_4();
			
			SimpleGapPenalty gapP = new SimpleGapPenalty();
			gapP.setOpenPenalty((short)5);
			gapP.setExtensionPenalty((short)2);
			
			psa = Alignments.getPairwiseAlignment(rec, ext, PairwiseSequenceAlignerType.GLOBAL, gapP, matrix);
			
			//System.out.println(psa);
			
			
		}
		catch(Exception e) {
			e.printStackTrace();
		}
		
		return psa;
		
//		List<AlignedSequence<DNASequence, NucleotideCompound>> aligned = null;
//		try {
//			//create a list of the DNASequences with the reconstructed sequences and the extant sequence
//			List<DNASequence> lst = new ArrayList<DNASequence>();
//			for(String x: recons){
//				lst.add(new DNASequence(x));
//			}
//			lst.add(new DNASequence(extant));
//			
//			
//			
////		//ybr seqs for test
////		//fastmlfreemarg
////		DNASequence sequence1 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCCATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCAATTGGTTTTGCAGACGGTATTTTTCAATTCCTTTTTTAGGTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCGTATCTTAA");
////		//fastML_free_joint
////		DNASequence sequence2 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCCATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTTGGTTTTGCAGACGGTATTTTTCAATTCCTTTTTTGGGTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCGTATCTTAA");
////		//fastmlspmarg
////		DNASequence sequence3 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCTATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTCGGTTTTGCAGACGGTATTTTTTAATTCCTTTTTTAGGTTTTGTTTCGTCTTTCCTTTTTTTTATTGTTTTCGTATCTTAA");
////		//fastML_sp_joint
////		DNASequence sequence4 = new DNASequence("GTGTCCCCTGTCTATATATATCCATTGACGGTCTATTCTATTCTCTTGGTACATGTTGAAGTGAGCGTTTTTTATTATTGCACTCGGTTTTGCAGACGGTATTTTTTAATTCCTTTTTTGGGTTTTGTTTCGTCTTTCCTTTTTTTTATTGTTTTCGTATCTTAA");
////		//prank_free
////		DNASequence sequence5 = new DNASequence("ATGTCCCATCTCTATATATATTCATTGACGG--CATTCTATTTTCTTGGT----ATTGAAATGAGCGTTTTTTATTATTACAATTGGTTTTGCAGACGGCACTTTCC--TTTCCTTTT---GTTTTGTTTCTTCTTTCCTTTTTTTTATTGTTCTCATATCTTAA");	
////		//pranksp
////		DNASequence sequence6 = new DNASequence("GAGGCTCCTGTCTATATATACCATTTGACATTCTATTCTATTCTACAACT----GTTAAAGTGAGTGGTTTTTATTGTAGAAGCTGGTTTTAAGGGTGGTATTTTTT--CTCCTTCCCTGGGTTTTGT------------ATATTTTATTGTTCTCCTATCTTAA");
////		//prequel_sp
////		DNASequence sequence7 = new DNASequence("GTGTTCCCTGTCTATATATATTCATTGACGGTATTTGTTTCTTTTGGTGTTGAAGTGAGCGTTTTTTATTATTGCACTTGGTTTTCCAGACGGTATTTTAATTTTATTTTTAAGTCTTGTTCTTCTTGCCTTTTTTTCTTTGTTCTGGTATCTTAA");	
////		
////		
////		List<DNASequence> lst = new ArrayList<DNASequence>();
////		lst.add(sequence1);
////		lst.add(sequence2);
////		lst.add(sequence3);
////		lst.add(sequence4);
////		lst.add(sequence5);
////		lst.add(sequence6);
////		lst.add(sequence7);
//		
//		Profile<DNASequence, NucleotideCompound> results = Alignments.getMultipleSequenceAlignment(lst);
//		
//		aligned = results.getAlignedSequences();
//		
//		/*
//		for(AlignedSequence<DNASequence, NucleotideCompound> x : huuuuh)
//		{
//			System.out.println(x.toString());
//		}
//		*/
//		}
//		catch(Exception e){
//			e.printStackTrace();
//		}
//		return aligned;
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
	
	//make a table with the reconstruction, the ORF, and the scores 
	public static void outTable(PrintWriter writer, int[] scores, int tool) throws IOException
	{
		//can take in an int array with every 2 steps corresponding to a new aligned/total codon fraction
		String[] ASRtools = new String[] {"FastML_free_marg", "FastML_free_joint", "FastML_sp_marg", "FastML_sp_joint", "prank_free", "prank_sp", "prequel_free", "prequel_sp"};	
		//if we did not find an ORF
		if(scores[0] == -1 && scores[1] == -1) {
			writer.printf("%18s | N/A \n", ASRtools[tool]);
		}
		else {
			writer.printf("%18s | %d/%d \n", ASRtools[tool], scores[0], scores[1]);			
		}
//		writer.println(locus + " aligned codons");
//		writer.println("____________________");
//		
//		//move through the array 
//		for(int i = 0; i < scores.length; i += 2)
//		{
//			
//			writer.format("%18s | %d/%d ",ASRtools[count], scores[i], scores[i+1]);		
//			count++;
//		}
		
		
		//writer.close();
		
	}
	
	public static String findORF(String seq, int startIndex)
	{
		int stopIndex = -1;
		String codon;
		//cycle through the string starting from the startIndex
		for(int i = startIndex; i < seq.length(); i+=3){
			if(i + 3 >= seq.length())
			{
				codon = seq.substring(i, seq.length());
			}
			else
			{
				codon = seq.substring(i, i+3);
			}
			if(isStop(codon))
			{
				stopIndex = i + 3;
				break;
			}
		}
		if(stopIndex == -1)
		{
			return "BAD";
		}
		else{
			return seq.substring(startIndex, stopIndex);
		}
	}
	
	public static ArrayList<Integer> findAllStart(String ORF)
	{
		//System.out.println("seq size: " + ORF.length());
        ArrayList<Integer> indexes = new ArrayList<Integer>();
        int index = 0;
        while(index != -1){
            index = ORF.indexOf(START, index + 3);
            if (index != -1) {
                indexes.add(index);
            }
        }
		/*TESTING
		System.out.println("Indexes array list size: " + indexes.size());
		System.out.println("seq size: " + ORF.length());
		for(int i = 0; i < indexes.size(); i++)
		{
			int j = indexes.get(i);

			System.out.println(j);
			System.out.println(ORF.substring(j,j+3));
			
		}
		*/
        return indexes;
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
    