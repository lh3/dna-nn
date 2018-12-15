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

//for (var i = 0; i < motif.length; ++i) print(motif[i]);

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
	if ((m = /^\s*\d+\s+\S+\s+\S+\s+\S+\s+(\S+)\s+(\d+)\s+(\d+)\s+\S+\s+[+C]\s+(\S+)\s+(\S+)/.exec(line)) != null) {
		var type = null;
		if (m[5] == "SINE/Alu") type = "1";
		else if (m[4] == "ALR/Alpha") type = "3";
		else if (m[4] == "BSR/Beta" || m[4] == "LSAU") type = "4";
		else if (m[4] == "HSATII") type = "2";
		else if ((m[5] == "Simple_repeat" || m[5] == "Satellite") && ((m4 = /^\(([ACGT]+)\)n/.exec(m[4])) != null)) {
			if (motif_hash[m4[1]] != null) {
				type = "2";
			} else if (m4[1].length % motif0.length == 0) {
				var c = 0;
				for (var j = 0; j < m4[1].length; j += motif0.length)
					if (motif_hash[m4[1].substr(j, j + motif0.length)] != null)
						++c;
				if (c * motif0.length == m4[1].length)
					type = "2";
			}
		}
		if (type != null) print(m[1], parseInt(m[2]) - 1, m[3], type, m[4], m[5]);
	}
}

buf.destroy();
file.close();
