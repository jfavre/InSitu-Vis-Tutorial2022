#include <mpi.h>
#include <adios2.h>
#include <thread>
#include <fides/DataSetReader.h>

#include <vtkm/cont/Algorithm.h>
#include <vtkm/cont/Initialize.h>
#include <vtkm/cont/ArrayRangeCompute.h>
#include <vtkm/worklet/DispatcherMapTopology.h>
#include <vtkm/worklet/WorkletMapTopology.h>
#include <vtkm/worklet/ScatterPermutation.h>
#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Scene.h>
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/View3D.h>
#include <vtkm/filter/contour/Contour.h>
#include <vtkm/filter/contour/ClipWithImplicitFunction.h>
#include <vtkm/cont/ColorTable.h>
#include <vtkm/io/VTKDataSetWriter.h>


using FieldInfoType = fides::metadata::Vector<fides::metadata::FieldInformation>;

int main(int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  auto opts = vtkm::cont::InitializeOptions::DefaultAnyDevice;
  vtkm::cont::InitializeResult config = vtkm::cont::Initialize(argc, argv, opts);
  adios2::ADIOS adios(MPI_COMM_WORLD);
  const std::string source_name = "source";

  fides::io::DataSetReader fides_reader("diffusion-catalyst-fides.json");
  std::unordered_map<std::string, std::string> paths;
  paths[source_name] = std::string("diffusion.bp");
  fides::DataSourceParams params;
  params["engine_type"] = "SST";
  fides_reader.SetDataSourceParameters(source_name, std::move(params));

  size_t step = 0;
  while(true) {
    auto status = fides_reader.PrepareNextStep(paths);
    fides::metadata::MetaData metaData = fides_reader.ReadMetaData(paths);
    if (status == fides::StepStatus::EndOfStream)
      {
      break;
      }

    fides::metadata::MetaData selections;
    fides::metadata::Vector<size_t> blockSelection;
    blockSelection.Data.push_back(0);

    FieldInfoType fieldSelection;
    fieldSelection.Data.push_back(fides::metadata::FieldInformation("temperature",
                                  vtkm::cont::Field::Association::Points));
    selections.Set(fides::keys::FIELDS(), fieldSelection);

    vtkm::cont::PartitionedDataSet output = fides_reader.ReadDataSet(paths, selections);
    auto inputData = output.GetPartition(0);
    if (!inputData.HasField("temperature", vtkm::cont::Field::Association::Points))
      {
      std::cerr << "Error: expected a temperature array. Did not get it." << std::endl;
      }
    else
      {
      const auto& scalarField = inputData.GetField("temperature");
      const auto& scalarHandle = scalarField.GetData().AsArrayHandle<vtkm::cont::ArrayHandle<double>>();
      vtkm::cont::ArrayHandle<vtkm::Range> rangeArray = vtkm::cont::ArrayRangeCompute(scalarHandle);
      auto rangePortal = rangeArray.ReadPortal();
      if (rangePortal.Get(0).Min != 0)
        {
        std::cerr << "Unexpected temperature min. Got " << rangePortal.Get(0).Min << std::endl;
        }
      if (rangePortal.Get(0).Max > 1.0)
        {
        std::cerr << "Unexpected temperature max range. Got " << rangePortal.Get(0).Max << std::endl;
        }
      }

    const std::vector<vtkm::Float64> isovalues = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};
    vtkm::Float64 isovalue = 0.240858;
// Create an isoline filter
    vtkm::filter::contour::Contour contour;
    contour.SetActiveField("temperature");
    contour.SetFieldsToPass({ "temperature" });
    contour.SetNumberOfIsoValues(7);
    //contour.SetIsoValues(isovalues);
    contour.SetIsoValue(0, 0.1);
    contour.SetIsoValue(1, 0.2);
    contour.SetIsoValue(2, 0.3);
    contour.SetIsoValue(3, 0.4);
    contour.SetIsoValue(4, 0.5);
    contour.SetIsoValue(5, 0.6);
    contour.SetIsoValue(6, 0.7);
        
    vtkm::cont::DataSet outputData = contour.Execute(inputData);
    std::string vtkfilename = "/dev/shm/foo." + std::to_string(step) + ".vtk";
    vtkm::io::VTKDataSetWriter writer(vtkfilename);
    writer.WriteDataSet(outputData);
    //outputData.PrintSummary(std::cout);
        
    vtkm::Bounds coordsBounds = inputData.GetCoordinateSystem().GetBounds();
    vtkm::rendering::Camera camera;
    camera.ResetToBounds(coordsBounds);

    camera.SetLookAt(vtkm::make_Vec(0.5f, 0.5f, 0.0f));
    camera.SetViewUp(vtkm::make_Vec(0.0f, 1.0f, 0.0f));
    //camera.SetClippingRange(1.f, 100.f);
    camera.SetFieldOfView(30.f);

    camera.SetPosition(vtkm::make_Vec(0.5f, 0.5f, 3.35f));

    vtkm::cont::ColorTable colorTable("inferno");

    // Create a mapper, canvas and view that will be used to render the scene
    vtkm::rendering::Scene scene;
    vtkm::rendering::MapperRayTracer mapper;
    vtkm::rendering::CanvasRayTracer canvas(1024,1024);
    vtkm::rendering::Color bg(0.2f, 0.2f, 0.2f, 1.0f);

    // Render an image of the 2D grid surface
    scene.AddActor(vtkm::rendering::Actor(inputData.GetCellSet(),
                                          inputData.GetCoordinateSystem(),
                                          inputData.GetField("temperature"),
                                          colorTable));
                                              
    // Render an image of the output isosurface
    scene.AddActor(vtkm::rendering::Actor(outputData.GetCellSet(),
                                          outputData.GetCoordinateSystem(),
                                          outputData.GetField("temperature"),
                                          colorTable));
    vtkm::rendering::View3D view(scene, mapper, canvas, camera, bg);
    //view.Initialize();
    view.Paint();
    std::string filename = "diffusion_step_" + std::to_string(step++) + ".png";
    view.SaveAs(filename);
    //std::cout << __FILE__ << ":" << __LINE__ << ": File " << filename << " written, going to next iteration of infinite loop." << std::endl;
    }

  MPI_Finalize();
}
