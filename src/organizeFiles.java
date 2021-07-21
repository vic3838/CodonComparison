import java.util.*;
import java.io.*;
public class organizeFiles {

	//this script is to organize all the output files into one directory
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		File targetDir = new File("D:\\College\\TECBio\\CodonComparison\\alignedCodonsPlot");
		File seekerPath = new File("D:\\College\\TECBio\\CodonComparison");
		String locus;				//locus, title of the directory
		String codingAlignment;		//file name of ORF alignment
		String rawAlignment;		//filename of raw reconstruction alignments
		String fileName;			//name of the file in the directory, used to copy to other directory
		String dest;				//destination of the file
		//array of the ORF files, YAL, YBR, etc
		File [] emergingORFs = seekerPath.listFiles(new MyFileNameFilter());
		for(File x : emergingORFs)
		{		
			String str = x.getAbsolutePath();				//path to the ORF folder
			System.out.println(str);
			if(x.isDirectory()) {		//if we have a directory
				locus = str.substring(str.lastIndexOf("\\") + 1);		//get the name of the ORF will be used in several ops
				codingAlignment = str + "\\" + locus + "_alignedCodonsByLongestCodingSequence.txt";		//set to file name of one file to be copied
				rawAlignment = str + "\\" + locus + "_alignedCodonsByRawReconstruction.txt";
				
				System.out.println(codingAlignment);
				System.out.println(rawAlignment);
				//we have the paths to the files, now we copy them to the directory
				//D:\College\TECBio\alignedCodonsPlot\......
				
				//first, we need to create the destination file address
				fileName = codingAlignment.substring(codingAlignment.lastIndexOf("\\") + 1);			
				dest = "D:\\College\\TECBio\\alignedCodonsPlot\\" + fileName;
				System.out.println(dest);
				copyFileUsingStream( new File(codingAlignment), new File(dest) );
				
				fileName = rawAlignment.substring(rawAlignment.lastIndexOf("\\") + 1);
				dest = "D:\\College\\TECBio\\alignedCodonsPlot\\" + fileName;
				System.out.println(dest);
				System.out.println();
				copyFileUsingStream( new File(rawAlignment), new File(dest) );
		
			}
			
		}
		
		
	}
	
	private static void copyFileUsingStream(File source, File dest) throws IOException {
	    InputStream is = null;
	    OutputStream os = null;
	    try {
	        is = new FileInputStream(source);
	        os = new FileOutputStream(dest);
	        byte[] buffer = new byte[1024];
	        int length;
	        while ((length = is.read(buffer)) > 0) {
	            os.write(buffer, 0, length);
	        }
	    } finally {
	        is.close();
	        os.close();
	    }
	}
	
	//iterate through all the ORF directories
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

}
