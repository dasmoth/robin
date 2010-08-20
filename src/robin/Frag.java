package robin;

public class Frag {
	public final String chr;
	public final int min;
	public final int max;
	public final Label label;
	// Also record strand for non-CpG analysis?
	
	Frag(String chr, int min, int max, Label label) {
		this.chr = chr; this.min = min; this.max = max;
		this.label = label;
	}

	public WindowAssociation[] windows;
	public int length;
}
