var file = arguments.length? new File(arguments[0]) : new File();
var buf = new Bytes();
var re = /(\d+)([MIDNSH])/g;

var n_err_bases = 0, n_err_reads = 0, tot_reads = 0, n_chimeric = 0, n_chimeric_reads = 0, n_unmapped = 0, n_perfect = 0, n_clipped = 0, tot_clip = 0;
var last = null, last_nm = 0, last_seg = 0, last_clip = 0;

function count(last_nm, last_seg, last_clip)
{
	++tot_reads;
	tot_clip += last_clip;
	if (last_nm == 0 && last_clip == 0 && last_seg == 1) ++n_perfect;
	if (last_nm > 0) ++n_err_reads, n_err_bases += last_nm;
	if (last_clip != 0) ++n_clipped;
	if (last_seg == 0) ++n_unmapped;
	else if (last_seg > 1) ++n_chimeric_reads, n_chimeric += last_seg - 1;
}

while (file.readline(buf) >= 0) {
	var line = buf.toString();
	if (line.charAt(0) == '@') continue;
	var m, t = line.split("\t");
	var flag = parseInt(t[1]), nm = 0;
	var seg = (flag&4)? 0 : 1;
	if ((m = /NM:i:(\d+)/.exec(line)) != null)
		nm = parseInt(m[1]);
	var name = t[0] + '/' + (flag>>6&3);
	if (name != last) {
		if (last) count(last_nm, last_seg, last_clip);
		last = name, last_nm = nm, last_seg = seg, last_clip = 0;
		while ((m = re.exec(t[5])) != null)
			if (m[2] == 'S' || m[2] == 'H')
				last_clip += parseInt(m[1]);
	} else last_seg += seg, last_nm += nm;
}
if (last) count(last_nm, last_seg, last_clip);

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
