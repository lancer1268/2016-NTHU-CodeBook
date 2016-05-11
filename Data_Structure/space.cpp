
double norm(void) {
    return hypot(x,y);
}

double dis(Point a,Point b)
{
    return (a-b).norm();
}
