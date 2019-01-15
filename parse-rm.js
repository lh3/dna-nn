#!/usr/bin/env k8

var motif0 = "GGAAT";
var comp_tbl = { 'A':'T', 'T':'A', 'C':'G', 'G':'C' };
var motif = [motif0], motif_alt = [], motif_hash = {};

// reverse complement
for (var i = 0; i < motif.length; ++i) {
	var x = motif[i], y = "";
	for (var j = x.length - 1; j >= 0; --j) {
		y += comp_tbl[x[j]];
	}
	motif_alt.push(y);
}
for (var i = 0; i < motif_alt.length; ++i)
	motif.push(motif_alt[i]);

// rotate
motif_alt = [];
for (var i = 0; i < motif.length; ++i) {
	var x = motif[i];
	for (var j = 1; j < x.length; ++j)
		motif_alt.push(x.substr(j) + x.substr(0, j));
}
for (var i = 0; i < motif_alt.length; ++i)
	motif.push(motif_alt[i]);

for (var i = 0; i < motif.length; ++i) motif_hash[motif[i]] = i;

// print usage
if (arguments.length == 0) {
	print("Usage: k8 parse-rm.js <rmsk.out>");
	exit(1);
}

// parse
var file = new File(arguments[0]);
var buf = new Bytes();

while (file.readline(buf) >= 0) {
	var m, m4, line = buf.toString();
	var ctg = null, st = -1, en = -1, rep = null, fam = null;
	if ((m = /^\s*\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+[+C]\s+(\S+)\s+(\S+)/.exec(line)) != null) {
		ctg = m[1], st = parseInt(m[2]) - 1, en = m[3], rep = m[4], fam = m[5];
	} else if ((m = /^\d+(\t\d+){4}\t(\S+)\t(\d+)\t(\d+)\t\S+\t[+-]\t(\S+)\t(\S+)\t(\S+)/.exec(line)) != null) {
		ctg = m[2], st = m[3], en = m[4], rep = m[5];
		fam = m[6] == m[7]? m[6] : m[6] + '/' + m[7];
	}
	if (ctg != null) {
		var type = null;
		if (fam == "SINE/Alu") type = "4";
		else if (rep == "ALR/Alpha") type = "2";
		else if (rep == "BSR/Beta" || rep == "LSAU") type = "3";
		else if (rep == "HSATII") type = "1";
		else if (fam == "LINE/L1") type = "5";
		else if ((fam == "Simple_repeat" || fam == "Satellite") && ((m4 = /^\(([ACGT]+)\)n/.exec(rep)) != null)) {
			if (motif_hash[m4[1]] != null) {
				type = "1";
			} else if (m4[1].length % motif0.length == 0) {
				var c = 0;
				for (var j = 0; j < m4[1].length; j += motif0.length)
					if (motif_hash[m4[1].substr(j, j + motif0.length)] != null)
						++c;
				if (c * motif0.length == m4[1].length)
					type = "1";
			}
		}
		if (type != null) print(ctg, st, en, type, rep, fam);
	}
}

buf.destroy();
file.close();
