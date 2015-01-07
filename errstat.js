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

var f1 = arguments.length? new File(arguments[0]) : new File();
var f2 = arguments.length >= 2? new File(arguments[1]) : null;
var buf = new Bytes();
var st;

var n_err_bases = 0, n_err_reads = 0, tot_reads = 0, n_chimeric = 0, n_chimeric_reads = 0, n_unmapped = 0, n_perfect = 0, n_clipped = 0, tot_clip = 0;
var n1 = 0, n2 = 0, na = 0;
var st1, st2, last1 = null, last2 = null;

while ((st1 = read1(f1, buf, last1)) != null) {
	++tot_reads;
	tot_clip += st1.cliplen;
	if (st1.nm == 0 && st1.cliplen == 0 && st1.n_segs == 1) ++n_perfect;
	if (st1.nm > 0) ++n_err_reads, n_err_bases += st1.nm;
	if (st1.cliplen != 0) ++n_clipped;
	if (st1.n_segs == 0) ++n_unmapped;
	else if (st1.n_segs > 1) ++n_chimeric_reads, n_chimeric += st1.n_segs - 1;
	if (f2) {
		st2 = read1(f2, buf, last2);
		if (st2 == null) throw Error("the 2nd file has fewer reads");
		//if (st1.name != st2.name) throw Error("different read names: "+st1.name+ " vs "+st2.name);
		if (st1.nm != st2.nm || st1.cliplen != st2.cliplen || st1.n_segs != st2.n_segs) {
			var t;
			if (st1.nm <= st2.nm && st1.cliplen <= st2.cliplen && st1.n_segs == 1) {
				t = "1", ++n1;
			} else if (st2.nm <= st1.nm && st2.cliplen <= st1.cliplen && st2.n_segs == 1) {
				t = "2", ++n2;
			} else t = "a", ++na;
			print(t, st1.name, st1.n_segs, st1.cliplen, st1.nm, st2.n_segs, st2.cliplen, st2.nm);
		}
		last2 = st2.next;
	}
	last1 = st1.next;
}

print("# reads:             " + tot_reads);
print("# perfect reads:     " + n_perfect);
print("# unmapped reads:    " + n_unmapped);
print("# chimeric reads:    " + n_chimeric_reads);
print("# chimeric events:   " + n_chimeric);
print("# reads w/ base err: " + n_err_reads);
print("# error bases:       " + n_err_bases);
print("# clipped reads:     " + n_clipped);
print("# clipped bases:     " + tot_clip);
if (f2) {
	print("# better reads:      " + n1);
	print("# worse reads:       " + n2);
	print("# ambiguous reads:   " + na);
}

buf.destroy();
if (f2) f2.close();
f1.close();
