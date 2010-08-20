package robin;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

public class MeDIPSet {
	private static WindowAssociation[] EMPTY_WA_ARRAY = new WindowAssociation[0];
	
	private List<Frag> frags = new ArrayList<Frag>();
	private List<Window> windows = new ArrayList<Window>();
	private Set<Label> labels = new HashSet<Label>();
	private String chr = null;
	
	private final int windowSize;
	private final int windowStep;
	
	public MeDIPSet(int windowSize, int windowStep)
	{
		this.windowSize = windowSize;
		this.windowStep = windowStep;
	}
	
	public void addFragment(String chr, int min, int max, Label label)
	{
		chr = chr.intern();
		if (this.chr == null) {
			this.chr = chr;
		} else if (this.chr != chr) {
			throw new IllegalArgumentException("Can only process one chr at a time");
		}
		labels.add(label);
		Frag f = new Frag(chr, min, max, label);
		frags.add(f);
		
		List<WindowAssociation> waL = new ArrayList<WindowAssociation>();
		int minWindow = min/windowStep;
		for (int w = minWindow; ; ++w) {
			Window ww = getWindow(w);
			if (ww.min > max) {
				break;
			}
			int lap = Math.min(max, ww.max) - Math.max(min, ww.min) + 1;
			if (lap > 0) {
				waL.add(new WindowAssociation(ww, (1.0 * lap) / (max - min + 1)));
			}
		}
		f.windows = waL.toArray(EMPTY_WA_ARRAY);
	}
	
	private Window getWindow(int w) {
		while (windows.size() <= w) {
			windows.add(null);
		}
		Window ww = windows.get(w);
		if (ww == null) {
			ww = new Window(chr, (w * windowStep) + 1, (w * windowStep) + windowSize, w);
			windows.set(w, ww);
		}
		return ww;
	}
	
}
