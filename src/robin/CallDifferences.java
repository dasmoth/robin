package robin;

import java.io.File;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

import org.bjv2.util.cli.App;
import org.bjv2.util.cli.Option;

@App(overview="Call differences between sets of MeDIP-seq data", generateStub=true)
public class CallDifferences {
	private static final Label FG_LABEL = new Label("FG");
	private static final Label BG_LABEL = new Label("BG");
	
	private String fg;
	private String bg;
	private int windowSize = 500;
	private int windowStep = 200;
	private String chr;
	private int min = -1;
	private int max = -1;
	
	
	
	@Option(help="Chromosome to analyse")
	public void setChr(String chr) {
		this.chr = chr;
	}


	public void setMax(int max) {
		this.max = max;
	}


	public void setMin(int min) {
		this.min = min;
	}


	public void setWindowSize(int windowSize) {
		this.windowSize = windowSize;
	}


	public void setWindowStep(int windowStep) {
		this.windowStep = windowStep;
	}


	@Option(help="File of 'background' (normal, control) reads")
	public void setBg(String bg) {
		this.bg = bg;
	}


	@Option(help="File of 'foreground' (mutant, treatment, affected, etc.) reads")
	public void setFg(String fg) {
		this.fg = fg;
	}


	/**
	 * @param args
	 */
	public void main(String[] args) {
		MeDIPSet ms = new MeDIPSet(windowSize, windowStep);
		loadReads(fg, ms, FG_LABEL);
		loadReads(bg, ms, BG_LABEL);
	}

	private void loadReads(String fn, MeDIPSet medip, Label label)
	{
		SAMFileReader r = new SAMFileReader(new File(fn), new File(fn + ".bai"));
		r.setValidationStringency(SAMFileReader.ValidationStringency.LENIENT);
		int len = r.getFileHeader().getSequence(chr).getSequenceLength();
		
		int qmin = min < 0 ? 500000 : min;
		int qmax = max < 0 ? len-500000 : max;
		
		CloseableIterator<SAMRecord> i = r.query(chr, qmin, qmax, true);
		while(i.hasNext()) {
			SAMRecord sr = i.next();
			if (sr.getProperPairFlag() && sr.getMappingQuality() >= 10) {
				int as = sr.getAlignmentStart();
				int ae = sr.getAlignmentEnd();
				int ms = sr.getMateAlignmentStart();
				
				if (as < ms) {
					// System.out.printf("%s\tSubsetBAM\tpair\t%d\t%d\t%d\t+\t.%n", sr.getReferenceName(), as, ms + (ae - as), sr.getMappingQuality());
					medip.addFragment(sr.getReferenceName(), as, ms + (ae - as), label);
				} else {
					// System.out.printf("%s\tSubsetBAM\tpair\t%d\t%d\t%d\t+\t.%n", sr.getReferenceName(), ms, ae, sr.getMappingQuality());
					medip.addFragment(sr.getReferenceName(), ms, ae, label);
				}
			}
		}
	}
}
