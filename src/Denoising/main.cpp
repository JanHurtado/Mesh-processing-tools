#include <QApplication>
#include <fstream>

#include <stdio.h>

#include "noise.h"
#include "iomesh.h"
#include "metric.h"
#include "color.h"
#include "denoising.h"

#include <ctime>

/*TriMesh uniFormLaplacian(TriMesh & mesh,vector<num_t> parameters)
{
    return uniformLaplacian(mesh,(int)(parameters[0]),0.5);
}

void parameterOptimization_UniformLaplacian(TriMesh & mesh)
{
    vector<num_t> parameters;
    vector<num_t> upperbounds;
    vector<num_t> lowerbounds;
    vector<num_t> deltas;
    vector<num_t> gradient;
    //num iterations
    parameters.push_back(1);
    upperbounds.push_back(50);
    lowerbounds.push_back(0);
    deltas.push_back(1);
    gradient.push_back(0);

    TriMesh current_mesh = mesh;
    TriMesh gradient_mesh = mesh;
    for(int i=0;i<max_iterations;i++)
    {
        current_mesh = uniformLaplacian(parameters);
        for(int j=0;j<parameters.size();j++)
        {
            num_t temp = parameters[j] + deltas[j];
            if(temp>=lowerbounds[j] && temp<=upperbounds[j])
            {
                vector<num_t> temp_parameters = parameters;
                temp_parameters[j]=temp;
                TriMesh temp_mesh = uniformLaplacian(mesh,temp_parameters);
                gradient[j] = temp;
            }
            else
                gradient[j] = 0;
        }
    }
}*/


/////////////////////////////////////

#include "iomesh.h"
#include "denoising.h"


/*int main(int argc, char *argv[])
{
    string inputFileName;
    string outputFileName;
    string algorithm_key_word ;
    vector<num_t> parameters;
    if(argc<2)
    {
        printf("No input file specified\n");
        printUsage();
        return 0;
    }
    else if(argc<3)
    {
        printf("No output file specified\n");
        printUsage();
        return 0;
    }
    else if(argc<4)
    {
        printf("No denoising algorithm specified\n");
        printKeyWords();
        printUsage();
        return 0;
    }
    inputFileName = argv[1];
    outputFileName = argv[2];
    algorithm_key_word = argv[3];
    for(int i=4;i<argc;i++)
    {
        num_t param;
        sscanf(argv[i], "%lf", &param);
        parameters.push_back(param);
    }
    TriMesh inputMesh,outputMesh;

    //inputFileName = "C:\\Users\\hurtado\\Documents\\code\\build\\src\\Debug\\noisy.off";
    if(!importMesh(inputMesh,inputFileName))
    {
        printf("Can't read input file \n");
        return 0;
    }

    if(algorithm_key_word == uniform_key_word)
    {
        AlgorithmUniformLaplacian algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == taubin_key_word)
    {
        AlgorithmTaubin algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == hc_key_word)
    {
        AlgorithmHCLaplacian algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == bilateral_key_word)
    {
        AlgorithmBilateral algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == fast_key_word)
    {
        AlgorithmFastAndEffective algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == bilateralnormal_key_word)
    {
        AlgorithmBilateralNormal algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
    }
    else if(algorithm_key_word == guided_key_word)
    {
        AlgorithmGuided algorithm(inputMesh);
        algorithm.setParameters(parameters);
        outputMesh = algorithm.getOutput();
        return 0;
    }
    else
    {
        printf("Not a valid algorithm \n");
        printKeyWords();
        return 0;
    }

    if(!exportMesh(outputMesh,outputFileName))
    {
        printf("Can't write output file \n");
        return 0;
    }
    else
        return 0;
}*/

void write1(ofstream & os)
{

}


void test()
{
    string input_dir = "C:/Users/hurtado/Documents/testCasesNew/";
    string output_dir = "C:/Users/hurtado/Documents/testCasesNewRes/";
    string input_file_extension = ".stl";
    string output_file_extension = ".stl";
    vector<string> fileNames = {"Testcase01Gerson","Testcase02Gerson","Testcase03Gerson","Testcase04Gerson","Testcase05Gerson","Testcase06Gerson", "Testcase07Gerson","Testcase08Gerson","Testcase09Gerson","Testcase10Gerson","Testcase11Gerson","Testcase12Gerson","Testcase13Gerson","Testcase14Gerson","Testcase15Gerson"};
    vector<DenoisingAlgorithm*> denoisingAlgorithms = {new AlgorithmUniformLaplacian(),new AlgorithmTaubin(),new AlgorithmHCLaplacian(),new AlgorithmBilateral(),new AlgorithmFastAndEffective(),new AlgorithmBilateralNormal(), new AlgorithmGuided()};
    vector<string> denoisingAlgorithmLabels = {uniform_key_word,taubin_key_word,hc_key_word,bilateral_key_word,fast_key_word,bilateralnormal_key_word,guided_key_word};
    vector<string> denoisingAlgorithmNames = {"Uniform",
                                             "Taubin",
                                             "HC Laplacian",
                                             "Bilateral",
                                             "Fast and Effective",
                                             "Bilateral Normal",
                                             "Guided"};
    vector<vector<num_t> > denoisingAlgorithmParameters = {uniform_parameters_default,taubin_parameters_default,hc_parameters_default,bilateral_parameters_default,fast_parameters_default,bilateralnormal_parameters_default,guided_parameters_default};
    vector<num_t (*)(TriMesh&,kd_tree_t&, kd_tree_t&,vector<num_t>&,TriMesh&,vector<num_t>&)> metrics = {&getMeanVertexDistanceError,&getL2VertexBasedError,&getMeanQuadricError,&getMeanSquareAngleError2,&getL2NormalBasedError,&getMeanTangentialError,&getMeanDiscreteCurvatureError,&getAreaError,&getVolError};
    vector<string> metricNames = {"Mean distance error","L2 vertex-based error metric","Quadric error metric","MSAE","L2 normal-based error metric","Tangential error metric","Discrete curvature error metric","Area error","Volume error"};
    for (int i=0;i<fileNames.size();i++)
    {
        TriMesh original_mesh;
        importMesh(original_mesh,input_dir+fileNames[i]+input_file_extension);
        TriMesh noisy_mesh = addNoise(original_mesh);
        exportMesh(noisy_mesh,output_dir+"noisy_"+fileNames[i]+output_file_extension);
        PointCloud v_cloud,f_cloud;
        generateVertexPointCloud(original_mesh,v_cloud);
        generateFaceCentroidPointCloud(original_mesh,f_cloud);
        kd_tree_t   index_vertices(3 /*dim*/, v_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_vertices.buildIndex();
        kd_tree_t   index_faces(3 /*dim*/, f_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_faces.buildIndex();
        vector<num_t> original_mesh_curvature;
        getVertexMeanCurvature(original_mesh,original_mesh_curvature);
        vector<vector<num_t> > res;
        ofstream fsT(output_dir+fileNames[i]+"_resT.txt");
        for(int j=0;j<denoisingAlgorithms.size();j++)
        {
            vector<num_t> temp_res;
            denoisingAlgorithms[j]->setInput(noisy_mesh);
            denoisingAlgorithms[j]->setParameters(denoisingAlgorithmParameters[j]);
            denoisingAlgorithms[j]->run();
            TriMesh denoised_mesh = denoisingAlgorithms[j]->getOutput();
            exportMesh(denoised_mesh,output_dir+denoisingAlgorithmLabels[j]+"_"+fileNames[i]+output_file_extension);
            vector<num_t> denoised_mesh_curvature;
            getVertexMeanCurvature(denoised_mesh,denoised_mesh_curvature);
            fsT<<denoisingAlgorithmLabels[j];
            for(int k=0;k<metrics.size();k++)
            {
                num_t val = (*(metrics[k]))(original_mesh,index_vertices,index_faces,original_mesh_curvature,denoised_mesh,denoised_mesh_curvature);
                temp_res.push_back(val);
                fsT<<" & "<<val;

            }
            fsT<<" \\\\"<<endl;
            res.push_back(temp_res);
        }
        fsT.close();
        ofstream fsC(output_dir+fileNames[i]+"_resC.txt");
        int num_metrics = metrics.size();
        int num_algorithms = denoisingAlgorithms.size();
        for(int ii=0;ii<num_metrics;ii++)
        {
            fsC<<endl<<metricNames[ii]<<endl<<endl;
            for(int jj=0;jj<num_algorithms;jj++)
            {
                num_t val = res[jj][ii];
                fsC<<"  ("<<val<<",{"<<denoisingAlgorithmNames[jj]<<"})"<<endl;
            }
            fsC<<endl<<endl;
        }
        fsC.close();
    }
}

void test2()
{
    string filename = "C:/Users/hurtado/Documents/testCasesNew/Testcase06Gerson.stl";
    string filenameOutput = "C:/Users/hurtado/Documents/testCasesNewRes/denoised.stl";
    string filenameNoisy = "C:/Users/hurtado/Documents/testCasesNewRes/noisy.stl";
    TriMesh mesh;
    importMesh(mesh,filename);
    TriMesh noisy = addNoise(mesh);
    TriMesh denoised;
    DenoisingAlgorithm * alg = new AlgorithmBilateral();
    alg->setInput(noisy);
    alg->setDefaultParameters();
    alg->run();

    //denoised = addNoise(mesh);
    denoised = alg->getOutput();
    exportMesh(denoised,filenameOutput);
    exportMesh(noisy,filenameNoisy);
}

int main(int argc, char *argv[])
{
    clock_t begin = clock();
    test();
    clock_t end = clock();
    double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
    printf("Time taken: %.2fs\n", elapsed_secs);
    return 0;
}



/*int mainn(int argc, char *argv[])
{
    string dir = "C:/Users/hurtado/Documents/testCases/";
    string filename = "testcase01.off";
    TriMesh mesh;
    importMesh(mesh,dir+filename);

    vector<num_t> parameters1 = {1};
    vector<num_t> parameters2 = {1,0.5,-0.53};
    vector<num_t> parameters3 = {1,0.5,0.5};
    vector<num_t> parameters4 = {20};
    //vector<num_t> parameters5 = {10,20,0.5};
    vector<num_t> parameters5 = {10,20};
    vector<num_t> parameters6 = {10,20,1,0.35};
    vector<vector<num_t> > parameters;
    parameters.push_back(parameters1);
    parameters.push_back(parameters2);
    parameters.push_back(parameters3);
    parameters.push_back(parameters4);
    parameters.push_back(parameters5);
    parameters.push_back(parameters6);
    vector<DenoisingAlgorithm*> algorithms;
    algorithms.push_back(new AlgorithmUniformLaplacian(mesh));
    algorithms.push_back(new AlgorithmTaubin(mesh));
    algorithms.push_back(new AlgorithmHCLaplacian(mesh));
    algorithms.push_back(new AlgorithmBilateral(mesh));
    algorithms.push_back(new AlgorithmFastAndEffective(mesh));
    algorithms.push_back(new AlgorithmBilateralNormal(mesh));
    for(int i=0;i<algorithms.size();i++)
    {
        algorithms[i]->setParameters(parameters[i]);
        algorithms[i]->run();
        TriMesh res = algorithms[i]->getOutput();
        exportMesh(res,dir+to_string(i)+filename);
    }
}*/

int mainnn(int argc, char *argv[])
{
    //QApplication::setStyle(QStyleFactory::create("cleanlooks"));
    //QApplication a(argc, argv);
   //MainWindow w;
    //w.show();
    //string dir = "C:/Users/hurtado/Documents/test/";
    string dir = "C:/Users/hurtado/Documents/testCases/";
    //string dir = "C:/Users/hurtado/Documents/test3/";
    //string dir = "C:/Users/hurtado/Documents/test2/";
    vector<string> filenames;
    //filenames.push_back("mesh.off");
    //filenames.push_back("mesh2.off");
    //filenames.push_back("mesh3.off");
    //filenames.push_back("mesh4.off");
    //filenames.push_back("mesh5.off");
    filenames.push_back("testcase01.off");
    /*filenames.push_back("testcase02.off");
    filenames.push_back("testcase03.off");
    filenames.push_back("testcase04.off");
    filenames.push_back("testcase05.off");
    filenames.push_back("testcase06.off");
    filenames.push_back("testcase08.off");
    filenames.push_back("testcase09.off");
    filenames.push_back("testcase10.off");
    filenames.push_back("testcase11.off");
    filenames.push_back("testcase12.off");
    filenames.push_back("testcase13.off");
    filenames.push_back("testcase14.off");
    filenames.push_back("testcase15.off");*/


    //filenames.push_back("block.off");
    //filenames.push_back("fandisk.off");
    //filenames.push_back("sharpSphere.off");
    //filenames.push_back("twelve.off");

    //filenames.push_back("bunny.off");
    //filenames.push_back("julius.off");
    //filenames.push_back("nicolo.off");
    ofstream res("res.txt");
    for (int i=0;i<filenames.size();i++)
    {
        TriMesh original_mesh;
        importMesh(original_mesh,dir+filenames[i]);
        TriMesh noisy_mesh = addNoise(original_mesh);
        //TriMesh denoised_uniform = uniformLaplacian(noisy_mesh, 10,1);
        TriMesh denoised_uniform = HildebrandtAndPolthier(noisy_mesh,2,0.5,10,1);
        //TriMesh denoised_uniform = original_mesh;
        cout<<"Filtros YA"<<endl;

        PointCloud v_cloud,f_cloud;
        generateVertexPointCloud(original_mesh,v_cloud);
        generateFaceCentroidPointCloud(original_mesh,f_cloud);
        cout<<"point clouds generadas YA"<<endl;
        kd_tree_t   index_vertices(3 /*dim*/, v_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_vertices.buildIndex();

        kd_tree_t   index_faces(3 /*dim*/, f_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_faces.buildIndex();
        cout<<"Indices construidos YA"<<endl;
        cout<<endl;
        cout<<getMeanVertexDistanceError(original_mesh,index_vertices,index_faces,denoised_uniform)<<endl;
        cout<<"Mean Vertex Distance Error  YA"<<endl;
        cout<<getL2VertexBasedError(original_mesh,index_vertices,index_faces,denoised_uniform)<<endl;
        cout<<"L2 Vertex Distance Error  YA"<<endl;
        cout<<getMeanQuadricError(original_mesh,index_vertices,denoised_uniform)<<endl;
        cout<<"Quadric Error  YA"<<endl;
        cout<<getMeanSquareAngleError2(original_mesh,index_vertices,index_faces,denoised_uniform)<<endl;
        cout<<"MSAE2 YA"<<endl;
        cout<<getMeanSquareAngleError(original_mesh,denoised_uniform)<<endl;
        cout<<"MSAE YA"<<endl;
        cout<<getL2NormalBasedError(original_mesh,index_vertices,index_faces,denoised_uniform)<<endl;
        cout<<"L2 normal error YA"<<endl;
        cout<<getMeanTangentialError(original_mesh,index_vertices,denoised_uniform)<<endl;
        cout<<"mean tangential error YA"<<endl;

        vector<num_t> mean;
        vector<num_t> gaussian;
        vector<num_t> principal;
        getAllCurvatures(original_mesh,mean,gaussian,principal);
        /*for(int j=0;j<mean.size();j++)
            cout<<mean[j]<<" "<<gaussian[j]<<" "<<principal[j]<<endl;*/
        //setPropertyColor(original_mesh,principal);
        exportMesh(original_mesh,dir+"original_"+filenames[i]);
        exportMesh(denoised_uniform,dir+"denoised_"+filenames[i]);

    }

    return 0;
}

int main2(int argc, char *argv[])
{
    string fileName = "C:/Users/hurtado/Documents/test4/chanvese3d_50Baby_BoederPC_Surface_NormalsAscii_SSD.ply";
    string fileName2 = "C:/Users/hurtado/Documents/test4/res.ply";
    string fileName3 = "C:/Users/hurtado/Documents/test4/chanvese3d_50Baby_BoederPC_Surface_NormalsAscii_SSD_Ascii2.ply";
    TriMesh mesh;
    importMesh(mesh,fileName3);
    for (TriMesh::HalfedgeIter h_it = mesh.halfedges_begin();h_it!=mesh.halfedges_end();h_it++)
    {
        //h_it->
    }
    vector<TriMesh::Normal> normals;
    getAllFaceNormals(mesh,normals);
    exportMesh(mesh,fileName2);

    return 0;
}

int main3(int argc, char *argv[])
{
    //QApplication::setStyle(QStyleFactory::create("cleanlooks"));
    //QApplication a(argc, argv);
   //MainWindow w;
    //w.show();
    //string dir = "C:/Users/hurtado/Documents/test/";
    string dir = "C:/Users/hurtado/Documents/testCases/";
    //string dir = "C:/Users/hurtado/Documents/test3/";
    //string dir = "C:/Users/hurtado/Documents/test2/";
    vector<string> filenames;
    //filenames.push_back("mesh.off");
    //filenames.push_back("mesh2.off");
    //filenames.push_back("mesh3.off");
    //filenames.push_back("mesh4.off");
    //filenames.push_back("mesh5.off");
    filenames.push_back("testcase01.off");
    /*filenames.push_back("testcase02.off");
    filenames.push_back("testcase03.off");
    filenames.push_back("testcase04.off");
    filenames.push_back("testcase05.off");
    filenames.push_back("testcase06.off");
    filenames.push_back("testcase08.off");
    filenames.push_back("testcase09.off");
    filenames.push_back("testcase10.off");
    filenames.push_back("testcase11.off");
    filenames.push_back("testcase12.off");
    filenames.push_back("testcase13.off");
    filenames.push_back("testcase14.off");
    filenames.push_back("testcase15.off");*/


    //filenames.push_back("block.off");
    //filenames.push_back("fandisk.off");
    //filenames.push_back("sharpSphere.off");
    //filenames.push_back("twelve.off");

    //filenames.push_back("bunny.off");
    //filenames.push_back("julius.off");
    //filenames.push_back("nicolo.off");
    ofstream res("res.txt");
    for (int i=0;i<filenames.size();i++)
    {
        TriMesh original_mesh;
        importMesh(original_mesh,dir+filenames[i]);
        //TriMesh noisy_mesh = original_mesh;
        TriMesh noisy_mesh = addNoise(original_mesh);
        //multiple_radius - multiple_sigma_c - normal_iteration_number - sigma_s - vertex_iteration_number;
        TriMesh denoised_guided = guided(noisy_mesh, 2.0, 1.0, 40, 0.35, 20);
        TriMesh denoised_bilateralNormal = bilateralNormal(noisy_mesh, 1, 40, 0.35, 20);
        //TriMesh denoised_uniform = HCLaplacian(noisy_mesh, 50,0.5,0.5);
        //TriMesh denoised_uniform = taubin(noisy_mesh, 5,0.5,-0.53);
        //TriMesh denoised_uniform = uniformLaplacian(noisy_mesh, 5,1);
        TriMesh denoised_uniform = HildebrandtAndPolthier(noisy_mesh,2,0.5,10,0.001);
        TriMesh denoised_bilateral = bilateral(noisy_mesh, 13);
        cout<<"Filtros YA"<<endl;
        num_t original_volume = getVolume(original_mesh);

        PointCloud v_cloud,f_cloud;
        generateVertexPointCloud(original_mesh,v_cloud);
        generateFaceCentroidPointCloud(original_mesh,f_cloud);
        cout<<"point clouds generadas YA"<<endl;
        kd_tree_t   index_vertices(3 /*dim*/, v_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_vertices.buildIndex();

        kd_tree_t   index_faces(3 /*dim*/, f_cloud, KDTreeSingleIndexAdaptorParams(10 /* max leaf */) );
        index_faces.buildIndex();
        cout<<"Indices construidos YA"<<endl;

        getMeanVertexDistanceError(original_mesh,index_vertices,index_faces,denoised_bilateral);
        cout<<"Mean Vertex Distance Error  YA"<<endl;
        getMeanSquareAngleError2(original_mesh,index_vertices,index_faces,denoised_bilateral);
        cout<<"MSAE2 YA"<<endl;

        num_t bilateralNormalError = getMeanSquareAngleError(denoised_bilateralNormal,original_mesh);
        num_t bilateralError = getMeanSquareAngleError(denoised_bilateral,original_mesh);
        num_t uniformError= getMeanSquareAngleError(denoised_uniform,original_mesh);
        num_t guidedError= getMeanSquareAngleError(denoised_guided,original_mesh);

        num_t bilateralNormalVolError = getVolError(original_volume,denoised_bilateralNormal);
        num_t bilateralVolError = getVolError(original_volume,denoised_bilateral);
        num_t uniformVolError = getVolError(original_volume,denoised_uniform);
        num_t guidedVolError = getVolError(original_volume,denoised_guided);

        cout<<bilateralError;

        //vector<num_t> curvature(original_mesh.n_vertices(),0.0);
        //getVertexGaussianCurvature(mesh, curvature);

        /*getVertexMeanCurvature(original_mesh, curvature);
        setPropertyColor(original_mesh,curvature);
        getVertexMeanCurvature(noisy_mesh, curvature);
        setPropertyColor(noisy_mesh,curvature);
        getVertexMeanCurvature(denoised_bilateralNormal, curvature);
        setPropertyColor(denoised_bilateralNormal,curvature);
        getVertexMeanCurvature(denoised_bilateral, curvature);
        setPropertyColor(denoised_bilateral,curvature);*/

        exportMesh(original_mesh,dir+"original_"+filenames[i]);
        exportMesh(noisy_mesh,dir+"noisy_"+filenames[i]);
        exportMesh(denoised_bilateralNormal,dir+"bilateralNormal_"+filenames[i]);
        exportMesh(denoised_bilateral,dir+"bilateral_"+filenames[i]);
        exportMesh(denoised_uniform,dir+"uniform_"+filenames[i]);
        exportMesh(denoised_guided,dir+"guided_"+filenames[i]);

        res<<"MSAE "<<filenames[i]<<" (uniform - bilateral - bilateral_normal - guided):  ";
        res<<uniformError<<"  "<<bilateralError<<"  "<<bilateralNormalError<<"  "<<guidedError<<"  "<<endl;
        res<<"Vol "<<filenames[i]<<" (uniform - bilateral - bilateral_normal - guided):  ";
        res<<uniformVolError<<"  "<<bilateralVolError<<"  "<<bilateralNormalVolError<<"  "<<guidedVolError<<"  "<<endl;
    }


    //importMesh(mesh,"C:/Users/hurtado/Documents/mesh4.off");


    //return a.exec();
    return 0;
}
