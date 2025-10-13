BEGIN {
    IFS="\t"
    OFS="\t"
}

{
    if (match($0, /dx:i:(-1|0|1)/, m) && match($0, /qs:f:([0-9]+.[0-9]+)/, n)) {
        if (n[1] >= 0) {
            n_seq[m[1]] += 1
            n_bases[m[1]] += length($10)
        }
        else {
            n_seq[2] += 1
            n_bases[2] += length($10)
        }
    }
}

END {
    print "low quality reads", n_seq[2]
    print "low quality reads bases", n_bases[2]
    print "simplex reads", n_seq[0]
    print "simplex bases", n_bases[0]
    print "duplex reads", n_seq[1]
    print "duplex bases", n_bases[1]
    print "simplex-duplex-offspring reads", n_seq[-1]
    print "simplex-duplex-offspring bases", n_bases[-1]
    print "duplex rate", 2*n_bases[1]/n_bases[0]
}