// Usage: k8 errstat.js [ec1.sam] [ec2.sam]

var re = /(\d+)([MIDNSH])/g;

function read1(f, b, t)
{
	var lines = [], m, name;

	if (t == null) {
		while (f.readline(b) >= 0) {
			if (b[0] == 64) continue;
			t = b.toString().split("\t");
			break;
		}
		if (t == null) return null; // end of file
		t[1] = parseInt(t[1]);
	}
	name = t[0] + '/' + (t[1]>>6&3);
	lines.push(t);
	t = null;

	while (f.readline(b) >= 0) {
		t = b.toString().split("\t");
		t[1] = parseInt(t[1]);
		var s = t[0] + '/' + (t[1]>>6&3);
		if (s != name) break;
		lines.push(t);
	}

	var st = { name:name, next:t, n_segs:0, nm:0, cliplen:0 };
	t = lines[0];
	if ((t[1]&4) == 0) {
		while ((m = re.exec(t[5])) != null)
			if (m[2] == 'S' || m[2] == 'H')
				st.cliplen += parseInt(m[1]);
	}
	for (var i = 0; i < lines.length; ++i) {
		t = lines[i];
		if (t[1]&4) continue;
		for (var j = 11; j < t.length; ++j)
			if (t[j].substr(0, 5) == "NM:i:")
				st.nm += parseInt(t[j].substr(5));
		++st.n_segs;
	}
	return st;
}

if (arguments.length < 2) {
	var file = arguments.length? new File(arguments[0]) : new File();
	var buf = new Bytes();
	var st;

	var n_err_bases = 0, n_err_reads = 0, tot_reads = 0, n_chimeric = 0, n_chimeric_reads = 0, n_unmapped = 0, n_perfect = 0, n_clipped = 0, tot_clip = 0;
	var last = null, st;

	while ((st = read1(file, buf, last)) != null) {
		++tot_reads;
		tot_clip += st.cliplen;
		if (st.nm == 0 && st.cliplen == 0 && st.n_segs == 1) ++n_perfect;
		if (st.nm > 0) ++n_err_reads, n_err_bases += st.nm;
		if (st.cliplen != 0) ++n_clipped;
		if (st.n_segs == 0) ++n_unmapped;
		else if (st.n_segs > 1) ++n_chimeric_reads, n_chimeric += st.n_segs - 1;
		last = st.next;
	}

	buf.destroy();
	file.close();

	print("# reads:             " + tot_reads);
	print("# perfect reads:     " + n_perfect);
	print("# unmapped reads:    " + n_unmapped);
	print("# chimeric reads:    " + n_chimeric_reads);
	print("# chimeric events:   " + n_chimeric);
	print("# reads w/ base err: " + n_err_reads);
	print("# error bases:       " + n_err_bases);
	print("# clipped reads:     " + n_clipped);
	print("# clipped bases:     " + tot_clip);
} else {
	var f1 = new File(arguments[0]);
	var f2 = new File(arguments[1]);
	var buf = new Bytes();

	buf.destroy();
	f2.close();
	f1.close();
}
