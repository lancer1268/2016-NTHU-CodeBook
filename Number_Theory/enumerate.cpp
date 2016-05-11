
void all_divdown(const LL &n) { // all n/x
    for(LL a=1;a<=n;a=n/(n/(a+1))) {
        dosomething;
    }
}

LL Catalan(const LL &n) { return (2n)! / n! / n! / (n+1); }
