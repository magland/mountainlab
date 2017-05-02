### Spike sorting via MountainLab/MountainSort involves the following steps

1. Install and configure MountainLab
2. Prepare the raw data
3. Prepare a sorting project
4. Run the sorting
5. View the results

### A note about processing

MountainLab processing is defined in layers. This tutorial shows the recommended way to do the processing. However you are by no means restricted to doing it this way -- for example, the need to launch a processing daemon may seem like overkill. All processing routines may be traced back to simple executables. An explanation of the processing layers and different ways to invoke spike sorting is found [[todo: processing layers]].

### 1. Install and configure MountainLab

This involves installing prerequisites, cloning the repository, compiling, and editing mountainlab.user.json. Installation instructions are [found here](doc/installation.md). You should test the installation by running the example before proceeding.

Once installed do the following:

> cd mountainlab

> cp mountainlab.default.json mountainlab.user.json

Note that the mountainlab.default.json file is part of the repository and should not be edited. The settings in the mountainlab.user.json file may be edited, and they override any settings in the default file. 

Your mountainsort.user.json file should now look something like:
```json
{
  "general":{
    "temporary_path":"/tmp",
    "max_cache_size_gb":40
  },
  "mountainprocess":{
    "max_num_simultaneous_processes":2,
    "processor_paths":["mountainprocess/processors","user/processors","packages"]
  },
  "prv":{
    "local_search_paths":["examples"],
    "servers":[
      {"name":"datalaboratory","passcode":"","host":"http://datalaboratory.org","port":8080,"path":"/prv"},
      {"name":"river","passcode":"","host":"http://river.simonsfoundation.org","port":8080,"path":"/prv"},
      {"name":"localhost","passcode":"","host":"http://localhost","port":8080,"path":"/prv"}
    ]
  },
  "kron":{
    "dataset_paths":["mountainsort/example_datasets","user/datasets"],
    "pipeline_paths":["mountainsort/pipelines","/user/pipelines"],
    "view_program_paths":["mountainsort/view_programs","/user/view_programs"]
  }
}
```

You should at least change the following configuration parameters:

general.temporary_path (default="/tmp"). This is a very important setting as it specifies the folder where all the temporary files and intermediate processing files are stored. Put it on a disk with a lot of space. Ideally it would be on an SSD drive but that's not necessary. More about the temporary directory elsewhere

prv.local_search_paths (default=["examples"]). Add the full path of the base directory where your raw data reside. The system will search recursively for the raw data files. More on that below. For example, set it to ["examples","/path/to/prvdata"].

The other settings are described in [[todo: prv_system and processing_layers]].

### 2. Prepare the raw data

The raw timeseries data must first be converted to .mda format. This is [described here](mda_format.md). Next, put it somewhere in the /path/to/prvdata as configured in mountainlab.user.json. That directory will be searched recursively. The principle (as will be described in more detail) is to separate the huge raw data files from the rest of the analysis procedure.

### 3. Prepare a sorting project

The directory structure for a sorting project is as follows:

```bash
project_name
  datasets.txt
  pipelines.txt
  curation.script
  datasets
    dataset1_name
      raw.mda.prv
      params.json
      geom.csv
    dataset2_name
      ...
    ```
```

These files are defined as follows:

datasets.txt is an index to the datasets, one line per dataset. For example:
ds1 datasets/dataset1_name
ds2 datasets/dataset2_name
The first column is an abreviated name and the second column is the relative path to the dataset folder.

pipelines.txt is an index to the processing pipelines. For example:
ms2 mountainsort_002.pipeline --curation=curation.script
ms_nf3 mountainsort_001.pipeline --curation=curation.script --num_features=3 --num_features2=3
The first column is an abreviated name. The second column is the name of the processing pipeline. The mountainsort_002.pipeline is distributed with MountainLab and is found in mountainlab/mountainsort/pipelines. That directory is searched by default because of the settings in mountainlab.user.json. Then each line contains processing parameters.

curation.script is optional and contains rules for rejecting or tagging clusters based on metrics computed as part of the sorting pipeline. You can find an example in mountainlab/examples/003_kron_mountainsort/curation.script. To include it in your sorting pipeline, you must 

1. Have a curation.script in the same location as the pipelines.txt and datasets.txt
2. Add the curation flag, --curation=curation.script to the pipelines.txt as seen above.

The curation script and the cluster and cluster-pair metrics are described in greater detail in [[todo: metrics and MountainView]].

raw.mda.prv must be created using the prv-create utility as follows:

> prv-create /path/to/raw/timeseries_data.mda /path/to/datasets/dataset1_name/raw.mda.prv

Take a look at the contents of raw.mda.prv. It is a JSON text file with information about the checksum of the raw data file. Thus it is a universal pointer to the raw data.  Even if you change the name or location of the raw data, the system will find it, assuming it is found somewhere in the prv search path (from mountainlab.user.json).

geom.csv is optional and contains the 2D or 3D locations of the electrodes. It is a comma-separated text file where each line represents an electrode channel, and the columns are the geometric coordinates. For example a tetrode might have the following geom.csv:

```bash
0,0
-1,1
1,1
0,2
```

This file is used in conjunction with the adjacency_radius sorting parameter and determines the local electrode neighborhoods. If adjacency_radius=0, or the geom.csv is not present, then there is only one electrode neighborhood containing all the channels. This file is also used by the viewer for display.

params.json contains sorting parameters that are specific to the dataset. At a minimum it should contain the sample rate in Hz as follows:

```json
{"samplerate":30000}
```

You can also specify whether to look for positive spike peaks (sign=1) , negative (sign=-1), or both (sign=0), as follows:

```json
{"samplerate":30000,"sign":1}
```

### 4. Run the sorting

Before running the sorting you should start the processing daemon:

> mp-daemon-start [name]

where [name] is the name of the daemon. Usually you would have only one of these running, so it could be your user name.

Then launch the sorting using:

> kron-run ms2 ds1

The output will go into the outputs/ms2–ds1 folder. In particular you will get a firings.mda file, which is [described here](doc/mda_format.md).

Further description of the daemon is found [[todo: processing_layers]].

### 5. View the results

To view the results use the following command:

> kron-view results ms2 ds1

Don't forget the word "results". This launches the MountainView GUI. Further description of the GUI can be found [[todo: MountainView]].
