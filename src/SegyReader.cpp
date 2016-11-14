#include "SegyReader.h"
#include <assert.h>

#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkStructuredPointsReader.h>
#include <vtkVolumeTextureMapper3D.h>
#include <vtkColorTransferFunction.h>
#include <vtkRenderWindow.h>
#include <vtkAxesActor.h>
#include <vtkImageShiftScale.h>
#include <vtkPolyDataMapper.h>
#include <vtkFloatArray.h>
#include <vtkPointData.h>
#include <vtkArrayData.h>
#include <vtkNew.h>
#include <vtkObjectFactory.h>
#include <set>
#include <vtkInformationVector.h>
#include <vtkInformation.h>
#include <vtkCellData.h>


SegyReader::~SegyReader()
{
    for(auto trace : data)
        delete trace;
}

bool SegyReader::LoadFromFile(string path) {
    ifstream in(path, ifstream::binary);
    if (!in)
    {
        cout << "File not found:" << path << endl;
        return false;
    }


    readHeader(in);

    int traceStartPos = 3600;// traces start after 3200 + 400 file header
    while (true)
    {
        Trace* pTrace = new Trace();
        if (!traceReader.readTrace(traceStartPos, in, formatCode, pTrace))
            break;
        data.push_back(pTrace);
    }

    in.close();
    return true;
}


bool SegyReader::readHeader(ifstream &in) {
    formatCode = IOUtil::Instance()->readShortInteger(binaryHeaderBytesPos.FormatCode, in);
    sampleCountPerTrace =  IOUtil::Instance()->readShortInteger(binaryHeaderBytesPos.NumSamplesPerTrace, in);
    return true;
}

#define PIXEL_TYPE float

bool SegyReader::ExportData3D(vtkImageData *imageData)
{

    set<int> crosslineNumbers, inlineNumbers;
    int minCrossLineNumber;
    int maxCrossLineNumber;
    int crosslineNum;

    for(auto trace : data)
    {
        crosslineNumbers.insert(trace->crosslineNumber);
        crosslineNum = trace->crosslineNumber;
        minCrossLineNumber = minCrossLineNumber < crosslineNum ? minCrossLineNumber : crosslineNum;
        maxCrossLineNumber = maxCrossLineNumber > crosslineNum ? maxCrossLineNumber : crosslineNum;
        inlineNumbers.insert(trace->inlineNumber);
    }

    int crossLineNumberStep = 1;
    int crosslineNumberCount = (maxCrossLineNumber - minCrossLineNumber) / crossLineNumberStep + 1;


    float min_data = INT_MAX;
    float max_data = INT_MIN;

    for(auto trace : data)
    {
        for(auto m : trace->data)
        {
            if(m < min_data)
                min_data = m;
            if(m > max_data )
                max_data = m;
        }
    }



    imageData->SetDimensions(sampleCountPerTrace, crosslineNumberCount, 1 );

    int type = VTK_FLOAT;
    imageData->SetScalarType(type, imageData->GetInformation());
    imageData->SetNumberOfScalarComponents(1, imageData->GetInformation());
    imageData->AllocateScalars(type, 1);
    float *ptr=(float *)imageData->GetScalarPointer();


    for (int k = 0; k < sampleCountPerTrace; k++) {
            for (int i = 0; i < crosslineNumberCount; i++) {

                *(ptr + i * sampleCountPerTrace + k) =
                        data[i]->data[k];
            }
  }


    return true;
}


bool SegyReader::GetImageData(vtkImageData *imageData)
{
  int crosslineNum;
  int minCrossLineNumber = INT_MAX;
  int maxCrossLineNumber = INT_MIN;

  for(auto trace : data)
  {
      crosslineNum = trace->crosslineNumber;
      if (crosslineNum == 0)
                  break;
      minCrossLineNumber = minCrossLineNumber < crosslineNum ? minCrossLineNumber : crosslineNum;
      maxCrossLineNumber = maxCrossLineNumber > crosslineNum ? maxCrossLineNumber : crosslineNum;
  }

  int crossLineNumberStep = 1;
  int crosslineNumberCount = (maxCrossLineNumber - minCrossLineNumber) / crossLineNumberStep + 1;

  int type = VTK_FLOAT;
  imageData->SetDimensions(sampleCountPerTrace, crosslineNumberCount, 1 );
  imageData->SetScalarType(type, imageData->GetInformation());
  imageData->SetNumberOfScalarComponents(1, imageData->GetInformation());
  imageData->AllocateScalars(type, 1);

  float *ptr=(float *)imageData->GetScalarPointer();


  for (int k = 0; k < sampleCountPerTrace; k++)
  {
        for (int i = 0; i < crosslineNumberCount; i++)
        {
          *(ptr + i * sampleCountPerTrace + k) = data[i]->data[k];
        }
  }


  return true;
}


bool SegyReader::AddScalars(vtkPolyData* polyData)
{

  vtkSmartPointer<vtkFloatArray> cellData =
    vtkSmartPointer<vtkFloatArray>::New();
  cellData->SetName("trace");
  cellData->SetNumberOfComponents(1);

  int crossLineCount = data.size();

  cellData->Allocate(crossLineCount * sampleCountPerTrace);

  float min_data = INT_MAX;
  float max_data = INT_MIN;

  for(auto trace : data)
  {
      for(auto m : trace->data)
      {
          if(m < min_data)
              min_data = m;
          if(m > max_data )
              max_data = m;
      }
  }

      for (int k = 0; k < sampleCountPerTrace; k++)
      {
        for (int i = 0; i < data.size(); i++)
        {
        cellData->InsertValue( i * sampleCountPerTrace + k, 256.0 * (data[i]->data[k] - min_data)/(max_data-min_data) );
        }
      }

    polyData->GetCellData()->SetScalars(cellData);
    polyData->GetCellData()->SetActiveScalars("trace");


    return true;
}


#include <vtkXMLPolyDataWriter.h>

bool SegyReader::ExportData2D(vtkPolyData * polyData)
{
    vtkSmartPointer<vtkPoints> points =
            vtkSmartPointer<vtkPoints>::New();

    vtkSmartPointer<vtkFloatArray> textureCoordinates =
            vtkSmartPointer<vtkFloatArray>::New();
    textureCoordinates->SetNumberOfComponents(2);
    textureCoordinates->SetName("TextureCoordinates");

    for(int k=0; k<sampleCountPerTrace; k++)
    {
        for(int i=0;i < data.size(); i++)
        {
            auto trace = data[i];
            float x = trace->xCoordinate / 10000.0;
            float y = trace->yCoordinate / 10000.0;

            float z = k * 100.0 / sampleCountPerTrace;
            points->InsertNextPoint(x, y, z);
            textureCoordinates->InsertNextTuple2(k * 1.0 / sampleCountPerTrace, i * 1.0 / data.size());
        }
    }

// Create a cell array to store the quad in
    vtkSmartPointer<vtkCellArray> quads =
            vtkSmartPointer<vtkCellArray>::New();

    for (int k = 1; k < sampleCountPerTrace; k++)
    {
        for(int i=1;i < data.size(); i++)
        {
            auto trace = data[i];

            vtkSmartPointer<vtkPolygon> polygon =
                    vtkSmartPointer<vtkPolygon>::New();
            polygon->GetPointIds()->SetNumberOfIds(4); //make a quad

            int id1 = k * data.size() + i;
            int id2 = (k - 1) * data.size() + i;
            int id3 = (k - 1) * data.size() + i - 1;
            int id4 = k * data.size() + i - 1;
            polygon->GetPointIds()->SetId(0, id1);
            polygon->GetPointIds()->SetId(1, id2);
            polygon->GetPointIds()->SetId(2, id3);
            polygon->GetPointIds()->SetId(3, id4);
            quads->InsertNextCell(polygon);
        }
    }

    // vtkSmartPointer<vtkTexture> texture =
    //        vtkSmartPointer<vtkTexture>::New();

    // vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
    // ExportData3D(imageData);
    // std::cerr << "imageData " << imageData << std::endl;
    // texture->SetInputDataObject(imageData);

    polyData->SetPoints(points);
    polyData->SetPolys(quads);
    this->AddScalars(polyData);
    // polyData->GetPointData()->SetTCoords(textureCoordinates);

    vtkNew<vtkXMLPolyDataWriter> writer;
    writer->SetInputData(polyData);
    writer->SetFileName("test.vtp");
    writer->Write();
    return polyData;
}





