#include "color.h"

/*void myColorMapping(num_t value,int &r,int &g,int &b)
{
    num_t a=(1-value)/0.25;	//invert and group
    num_t X = int(a);
    num_t Y=int(255*(a-X));
    //var X=Math.floor(a);	//this is the integer part
    //var Y=Math.floor(255*(a-X)); //fractional part from 0 to 255
    switch(int(X))
    {
        case 0: r=255;g=Y;b=0;break;
        case 1: r=255-Y;g=255;b=0;break;
        case 2: r=0;g=255;b=Y;break;
        case 3: r=0;g=255-Y;b=255;break;
        case 4: r=0;g=0;b=255;break;
    }
}

void setPropertyColor(TriMesh &mesh, vector<num_t> & property)
{
    num_t max_value = *max_element(property.begin(),property.end());
    num_t min_value = *min_element(property.begin(),property.end());
    min_value = 0;
    max_value = 1;
    for (TriMesh::VertexIter v_it=mesh.vertices_begin(); v_it!=mesh.vertices_end(); ++v_it)
    {
        num_t value = property[v_it.handle().idx()];
        //cout<<value<<endl;
        //value = (abs(value));

        value = (value - min_value) / (max_value - min_value);
        num_t final_value = max(0.0,min(value,1.0));
        int r,g,b;
        myColorMapping(final_value,r,g,b);
        if (value > 1) value = 1;
        if (value < 0) value = 0;
        //mesh.set_color(v_it.handle(),TriMesh::Color(r,g,b));
        mesh.set_color(v_it.handle(),TriMesh::Color(final_value*255,final_value*255,final_value*255));
        //cout<<value<<endl;
    }
    //cout<<*max_element(property.begin(),property.end())<<endl;
    //cout<<*min_element(property.begin(),property.end())<<endl;
}*/
