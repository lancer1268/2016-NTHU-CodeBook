Circle outcircle(Point a,Point b,Point c)
{
    if(TwoPointCircle(a,b).incircle(c)) return TwoPointCircle(a,b);
    if(TwoPointCircle(b,c).incircle(a)) return TwoPointCircle(b,c);
    if(TwoPointCircle(c,a).incircle(b)) return TwoPointCircle(c,a);

    Point ret;
    double a1=b.x-a.x, b1=b.y-a.y, c1=(a1*a1+b1*b1)/2;
    double a2=c.x-a.x, b2=c.y-a.y, c2=(a2*a2+b2*b2)/2;
    double d = a1*b2 - a2*b1;
    ret.x = a.x + (c1*b2-c2*b1)/d;
    ret.y = a.y + (a1*c2-a2*c1)/d;
    return Circle(ret,Dis(ret,a));
}

Circle SmallestCircle(vector<Point> &p) // save Points in p
{
    int n = p.size();
    if(n==1) return Circle(p[0],0.0);
    if(n==2) return TwoPointCircle(p[0],p[1]);

    random_shuffle(p.begin(),p.end());
    Circle c = Circle(p[0]);
    for(int i=0;i<n;++i){
        if(c.incircle(p[i])) continue;
        c = Circle(p[i]);
        for(int j=0;j<i;++j){
            if(c.incircle(p[j])) continue;
            c = TwoPointCircle(p[i],p[j]);
            for(int k=0;k<j;++k){
                if(c.incircle(p[k])) continue;
                c = outcircle(p[i],p[j],p[k]);
            }
        }
    }
    return c;
}
