#!/usr/bin/env nodejs

var common=require('./common.js');
var basic_sort=require('./ms2_001.js');
var os=require('os');

exports.spec=function() {
	var spec0=basic_sort.spec();
	spec0.name='mountainsort.ms2_001_multichannel';
	spec0.version=spec0.version+'-0.15';
	spec0.inputs.push({name:'geom'});
	spec0.parameters.push({name:'adjacency_radius'});
	return common.clone(spec0);
};

exports.run=function(opts,callback) {
	var tmpfiles=[];
	opts.temp_prefix=opts.temp_prefix||'00';
	if (!opts._tempdir) {
		console.error('opts._tempdir is empty');
		process.exit(-1);
	}
	if (!opts.samplerate) {
		console.error('opts.samplerate is zero or empty');
		process.exit(-1);	
	}
	if (!opts.adjacency_radius) {
		console.error('opts.adjacency_radius is zero or empty');
		process.exit(-1);		
	}

	var info={};
	read_info_from_input_files(function() {
		var num_threads=Number(opts.num_threads)||os.cpus().length;
		var num_threads_within_neighborhood=Math.floor(num_threads/info.M);
		if (num_threads_within_neighborhood<1) num_threads_within_neighborhood=1;
		var num_parallel_neighborhoods=info.M;
		if (num_parallel_neighborhoods>num_threads) num_parallel_neighborhoods=num_threads;

		var nbhd_steps=[];
		var all_nbhd_firings=[];
		for (var ch=1; ch<=info.M; ch++) {
			add_neighborhood_step(ch);
		}
		function add_neighborhood_step(ch) {
			nbhd_steps.push(function(cb) {
				var channels=get_neighborhood_channels(info.geom_coords,ch,opts.adjacency_radius);
				var firings0=mktmp('firings_nbhd'+ch+'.mda');
				all_nbhd_firings.push(firings0);
				{
					var opts2=common.clone(opts);
					opts2.central_channel=ch;
					opts2.firings_out=firings0;
					opts2.cluster_metrics_out=null; //we'll compute this at the end
					opts2.num_threads=num_threads_within_neighborhood;
					opts2.consolidate_clusters='true';
					opts2.channels=channels.join(',');
					opts2._temp_prefix=(opts._temp_prefix||'00')+'-nbhd'+ch;
					console.log ('Running basic sort for electrode '+ch+'...\n');
					basic_sort.run(opts2,function() {
						cb();
					});
				}
			});
		}
		console.log ('Num neighborhood steps: '+nbhd_steps.length+', num_parallel_neighborhoods: '+num_parallel_neighborhoods);
		//Run all steps
		common.foreach(nbhd_steps,{num_parallel:num_parallel_neighborhoods},function(ii,step,cb) {
			console.log ('');
			console.log ('--------------------------- NEIGHBORHOOD '+(ii+1)+' of '+nbhd_steps.length +' -----------');
			var timer=new Date();
			step(function() {
				console.log ('Elapsed time for neighborhood '+(ii+1)+' (sec): '+get_elapsed_sec(timer));
				cb();
			});
		},function() {
			combine_firings();
		});

		function combine_firings() {
			var firings_combined=mktmp('firings_combined.mda');
			common.mp_exec_process('mountainsort.combine_firings',
					{firings_list:all_nbhd_firings},
					{firings_out:firings_combined},
					{increment_labels:'true'},
					function() {
						common.copy_file(firings_combined,opts.firings_out,function() {
							cleanup(function() {
							});
						});
					}
			);
		}
	});
	
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////
	////////////////////////////////////////////////////////


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ //

	function cleanup(callback) {
		common.remove_temporary_files(tmpfiles,callback);
	}

	function fit_stage(timeseries,firings,firings_out,callback) {
		common.mp_exec_process('mountainsort.fit_stage',
			{timeseries:timeseries,firings:firings},
			{firings_out:firings_out},
			{
			},
			callback
		);
	}

	function mktmp(name) {
		var temp_prefix=opts._temp_prefix||'00';
		var path=opts._tempdir+'/'+temp_prefix+'-'+name;
		tmpfiles.push(path);
		return path;
	}

	function read_info_from_input_files(callback) {
		console.log ('Reading info from input files...');
		info.geom_coords=read_geom_coords(opts.geom);
		if (typeof(opts.timeseries)!='object') {
			common.read_mda_header(opts.timeseries,function (header) { // Read the .mda header for the timeseries
				info.M=Number(header.dims[0]);
				info.N=Number(header.dims[1]);
				callback();
			});
		}
		else {
			get_header_for_concatenation_of_timeseries(opts.timeseries,function(header) {
				info.M=Number(header.dims[0]);
				info.N=Number(header.dims[1]);
				callback();
			});
		}
	}

	function get_header_for_concatenation_of_timeseries(ts_list,callback) {
		if (ts_list.length===0) {
			console.error('ts_list is empty.');
			process.exit(-1);
		}
		common.read_mda_header(ts_list[0],function(header0) {
			header0.dims[1]=0;
			common.foreach(ts_list,{},function(ii,ts,cb) {
				common.read_mda_header(ts,function(header1) {
					header0.dims[1]+=header1.dims[1];
					cb();
				});
			},function() {
				callback(header0);
			});
		});
		
	}
};

function get_elapsed_sec(timer) {
	var stamp=new Date();
	return (stamp-timer)/1000;
}

function get_neighborhood_channels(coords,chan,adjacency_radius) {
	var pt0=coords[chan-1];
	var ret=[];
	for (var i=0; i<coords.length; i++) {
		var dist0=compute_distance(coords[i],pt0);
		if (dist0<=adjacency_radius)
			ret.push(i+1);
	}
	return ret;
}

function compute_distance(pt1,pt2) {
	var sumsqr=0;
	for (var i=0; i<pt1.length; i++) {
		sumsqr=sumsqr+(pt1[i]-pt2[i])*(pt1[i]-pt2[i]);
	}
	return Math.sqrt(sumsqr);
}

function read_geom_coords(geom) {
	var ret=[];
	var txt=common.read_text_file(geom);
	var lines=txt.split('\n');
	for (var i=0; i<lines.length; i++) {
		var line=lines[i].trim();
		if (line) {
			var vals=lines[i].split(',');
			if (vals.length>0) {
				var pt=[];
				for (var j=0; j<vals.length; j++) {
					pt.push(Number(vals[j].trim()));
				}
				ret.push(pt);
			}
		}
	}
	return ret;
}