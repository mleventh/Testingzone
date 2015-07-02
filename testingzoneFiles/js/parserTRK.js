
goog.require('THREE')
goog.require('goog.math.Vec3');



/**
 * Create a parser for the binary .TRK format.
 *
 * @constructor
 */
parserTRK = function() {

  //
  // call the standard constructor of X.base
  //goog.base(this);

  //
  // class attributes

  /**
   * @inheritDoc
   * @const
   */
  this._classname = 'parserTRK';

};
// inherit from X.parser
//goog.inherits(X.parserTRK, X.parser);

console.log("noise")

parserTRK.prototype.parse = function(container, object, data, flag) {

  console.log("!!!!")
  //X.TIMER(this._classname + '.parse');

  this._data = data;
  console.log(data)

  //
  // call the standard constructor of X.base
  //goog.base(this);

  //
  // class attributes

  /**
   * @inheritDoc
   * @const
   */
  this._classname = 'parser';

  /**
   * The data.
   *
   * @type {?ArrayBuffer}
   * @protected
   */
  //this._data = null;

  /**
   * The pointer to the current byte.
   *
   * @type {!number}
   * @protected
   */
  this._dataPointer = 0;

  /**
   * The native endianness flag. Based on
   * https://github.com/kig/DataStream.js/blob/master/DataStream.js
   *
   * @type {!boolean}
   * @protected
   */
  this._nativeLittleEndian = new Int8Array(new Int16Array([ 1 ]).buffer)[0] > 0;

  /**
   * The data-specific endianness flag.
   *
   * @type {!boolean}
   * @protected
   */
  this._littleEndian = true;

  /**
   * The min value of the last parsing attempt.
   *
   * @type {!number}
   * @protected
   */
  this._lastMin = -Infinity;

  /**
   * The max value of the last parsing attempt.
   *
   * @type {!number}
   * @protected
   */
  this._lastMax = Infinity;



/**
 * Parse data and configure the given object. When complete, a
 * X.parser.ModifiedEvent is fired.
 *
 * @param {!X.base}
 *          container A container which holds the loaded data. This can be an
 *          X.object as well.
 * @param {!X.dataContainer}
 *          object The object to configure.
 * @param {!ArrayBuffer}
 *          data The data to parse.
 * @param {*}
 *          flag An additional flag.
 * @throws {Error}
 *           An exception if something goes wrong.
 */
parserTRK.prototype.parse = function(container, object, data, flag) {

  throw new Error('The function parse() should be overloaded.');
  
};

//
// PARSE FUNCTIONS
//
//
/**
 * Get the min and max values of an array.
 *
 * @param {Array|Uint8Array|Uint16Array|Uint32Array|null}
 *          data The data array to analyze.
 * @return {!Array} An array with length 2 containing the [min, max] values.
 */
parserTRK.prototype.arrayMinMax = function(data) {

  var _min = Infinity;
  var _max = -Infinity;

  // buffer the length
  var _datasize = data.length;

  var i = 0;
  for (i = 0; i < _datasize; i++) {

    if(!isNaN(data[i])) {

      var _value = data[i];
      _min = Math.min(_min, _value);
      _max = Math.max(_max, _value);

    }

  }

  return [ _min, _max ];
 
};

/**
 * Create a string from a bunch of UChars. This replaces a
 * String.fromCharCode.apply call and therefor supports more platforms (like the
 * Android stock browser).
 *
 * @param {!Array|Uint8Array}
 *          array The Uint8Array.
 * @param {?number=}
 *          start The start position. If undefined, use the whole array.
 * @param {?number=}
 *          end The end position. If undefined, use the whole array.
 * @return {string} The created string.
 */
parserTRK.prototype.parseChars = function(array, start, end) {

  // without borders, use the whole array
  if (start === undefined) {

    start = 0;

  }
  if (end === undefined) {

    end = array.length;

  }

  var _output = '';
  // create and append the chars
  var i = 0;
  for (i = start; i < end; ++i) {

    _output += String.fromCharCode(array[i]);

  }

  return _output;

};


/**
 * Jump to a position in the byte stream.
 *
 * @param {!number}
 *          position The new offset.
 */
parserTRK.prototype.jumpTo = function(position) {

  this._dataPointer = position;

};

parserTRK.scan = function(type, chunks) {

    if (!goog.isDefAndNotNull(chunks)) {

      chunks = 1;

    }

    var _chunkSize = 1;
    var _array_type = Uint8Array;

    switch (type) {

    // 1 byte data types
    case 'uchar':
      break;
    case 'schar':
      _array_type = Int8Array;
      break;
    // 2 byte data types
    case 'ushort':
      _array_type = Uint16Array;
      _chunkSize = 2;
      break;
    case 'sshort':
      _array_type = Int16Array;
      _chunkSize = 2;
      break;
    // 4 byte data types
    case 'uint':
      _array_type = Uint32Array;
      _chunkSize = 4;
      break;
    case 'sint':
      _array_type = Int32Array;
      _chunkSize = 4;
      break;
    case 'float':
      _array_type = Float32Array;
      _chunkSize = 4;
      break;
    case 'complex':
      _array_type = Float64Array;
      _chunkSize = 8;
      break;
    case 'double':
      _array_type = Float64Array;
      _chunkSize = 8;
      break;

    }
    
    // increase the data pointer in-place
    var _bytes = new _array_type(data.slice(this._dataPointer,
        this._dataPointer += chunks * _chunkSize));
        console.log(this._dataPointer)
    // if required, flip the endianness of the bytes
    if (this._nativeLittleEndian != this._littleEndian) {

      // we need to flip here since the format doesn't match the native endianness
      _bytes = this.flipEndianness(_bytes, _chunkSize);

    }

    if (chunks == 1) {

      // if only one chunk was requested, just return one value
      return _bytes[0];

    }

    // return the byte array
    return _bytes;
};

/**
 * Flips typed array endianness in-place. Based on
 * https://github.com/kig/DataStream.js/blob/master/DataStream.js.
 *
 * @param {!Object}
 *          array Typed array to flip.
 * @param {!number}
 *          chunkSize The size of each element.
 * @return {!Object} The converted typed array.
 */
parserTRK.prototype.flipEndianness = function(array, chunkSize) {

  var u8 = new Uint8Array(array.buffer, array.byteOffset, array.byteLength);
  for ( var i = 0; i < array.byteLength; i += chunkSize) {

    for ( var j = i + chunkSize - 1, k = i; j > k; j--, k++) {

      var tmp = u8[k];
      u8[k] = u8[j];
      u8[j] = tmp;

    }

  }

  return array;

};

/**
 * Compute RAS bonding box fron IJK dimensions.
 *
 * @param {!Float32Array} IJKToRAS The IJK to RAS transformation.
 * @param {!Array} MRIdim The IJK dimensions.
 *
 * @return The RAS bounding box.
 * @static
 */
parserTRK.computeRASBBox = function(IJKToRAS, MRIdim){

  var _rasBB = [Number.MAX_VALUE, -Number.MAX_VALUE,
               Number.MAX_VALUE, -Number.MAX_VALUE,
               Number.MAX_VALUE, -Number.MAX_VALUE];

  var ijkTarget = goog.vec.Vec4.createFloat32FromValues(0, 0, 0, 1);
  var rasResult = goog.vec.Vec4.createFloat32();
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(0, 0, MRIdim[2]-1, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(0, MRIdim[1]-1, 0, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(MRIdim[0]-1, 0, 0, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(MRIdim[0]-1, MRIdim[1]-1, 0, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(MRIdim[0]-1, 0, MRIdim[2]-1, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(0, MRIdim[1]-1, MRIdim[2]-1, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

  ijkTarget = goog.vec.Vec4.createFloat32FromValues(MRIdim[0]-1, MRIdim[1]-1, MRIdim[2]-1, 1);
  goog.vec.Mat4.multVec4(IJKToRAS, ijkTarget, rasResult);

  _rasBB[0] = rasResult[0] < _rasBB[0] ? rasResult[0] : _rasBB[0];
  _rasBB[1] = rasResult[0] > _rasBB[1] ? rasResult[0] : _rasBB[1];
  _rasBB[2] = rasResult[1] < _rasBB[2] ? rasResult[1] : _rasBB[2];
  _rasBB[3] = rasResult[1] > _rasBB[3] ? rasResult[1] : _rasBB[3];
  _rasBB[4] = rasResult[2] < _rasBB[4] ? rasResult[2] : _rasBB[4];
  _rasBB[5] = rasResult[2] > _rasBB[5] ? rasResult[2] : _rasBB[5];

return _rasBB;
}

/**
 * Create the IJK volume.
 *
 * @param {!Float32Array} _data The target Bounding Box.
 * @param {!Array} _dims The line origin.
 * @param {!number} _max The maximum intensity value.
 *
 * @return The IJK volume and the IJK 'normalized' volume.
 * @static
 */
parserTRK.createIJKVolume = function(_data, _dims, _max, _min){
  // initiate variables
  // allocate images
  
  var _image = new Array(_dims[2]);
  var _imageN = new Array(_dims[2]);
  var _nb_pix_per_slice = _dims[1] * _dims[0];
  var _pix_value = 0;
  var _i = 0;
  var _j = 0;
  var _k = 0;
  var _data_pointer = 0;

  for (_k = 0; _k < _dims[2]; _k++) {

    // get current slice
    var _current_k = _data.subarray(_k * (_nb_pix_per_slice), (_k + 1)
        * _nb_pix_per_slice);
    // initiate data pointer
    _data_pointer = 0;

    // allocate images
    _imageN[_k] = new Array(_dims[1]);
    _image[_k] = new Array(_dims[1]);

    for (_j = 0; _j < _dims[1]; _j++) {

      // allocate images
      _imageN[_k][_j] = new _data.constructor(_dims[0]);
      _image[_k][_j] = new _data.constructor(_dims[0]);
      for (_i = 0; _i < _dims[0]; _i++) {
        _pix_value = _current_k[_data_pointer];
        _imageN[_k][_j][_i] = 255 * ((_pix_value - _min) / (_max - _min));
        _image[_k][_j][_i] = _pix_value;
        _data_pointer++;

      }
    }
  }

  return [_image, _imageN];
};

/**
 * Compute intersection between line and a bounding box
 *
 * @param {!Array} _bbox The target Bounding Box.
 * @param {!Float32Array} _sliceOrigin The line origin.
 * @param {!Float32Array} _sliceNormal The line normal.
 *
 * @return The intersection points and the 'non-intersection' points.
 * @static
 */
parserTRK.intersectionBBoxLine = function(_bbox, _sliceOrigin, _sliceNormal){

  var _solutionsIn = new Array();
  var _solutionsOut = new Array();

  // xmin, xmax, ymin, ymax, zmin, zmax
  for(var _i = 0; _i < 6; _i++) {

    var _i2 = Math.floor(_i/2);
    var _i3 = (_i2 + 1)%3;
    var _i4 = (_i2 + 2)%3;
    var _j1 = (2 + (2*_i2))%6;
    var _j2 = (4 + (2*_i2))%6;
    var _dir = _i2;


    var _sol0 = _bbox[_i];
    var _invN1 = 1/_sliceNormal[_i2];

    var _t = (_sol0 - _sliceOrigin[_i2])*_invN1;

    // if _t infinity, we are //
    if(_t != Infinity && _t != -Infinity) {

      var _sol1 = _sliceOrigin[_i3] + _sliceNormal[_i3]*_t;
      var _sol2 = _sliceOrigin[_i4] + _sliceNormal[_i4]*_t;

      // in range?
      if( (_sol1 >= _bbox[_j1] && _sol1 <= _bbox[_j1+1]) &&
          (_sol2 >= _bbox[_j2] && _sol2 <= _bbox[_j2+1])) {

        var _sol = new Array();
        _sol[_i2] = _bbox[_i];
        _sol[_i3] = _sol1;
        _sol[_i4] = _sol2;

        _solutionsIn.push(_sol);

      }
      else {

        var _sol = new Array();
        _sol[_i2] = _bbox[_i];
        _sol[_i3] = _sol1;
        _sol[_i4] = _sol2;

        _solutionsOut.push(_sol);

      }
    }
  }

  return [_solutionsIn, _solutionsOut];
};

/**
 * Compute intersection between plane and a bounding box
 *
 * @param {!Array} _bbox The target Bounding Box.
 * @param {!Float32Array} _sliceOrigin The plane origin.
 * @param {!Float32Array} _sliceNormal The plane normal.
 *
 * @return The intersection points and the 'non-intersection' points.
 * @static
 */
parserTRK.intersectionBBoxPlane = function(_bbox, _sliceOrigin, _sliceNormal){

  var _solutionsIn = new Array();
  var _solutionsOut = new Array();

  // xmin, xmax, ymin, ymax, zmin, zmax
  for(var _i = 0; _i < 6; _i++) {
    //
    var _i2 = Math.floor(_i/2);
    var _i3 = (_i2 + 1)%3;
    var _i4 = (_i2 + 2)%3;
    var _j3 = (4 + (2*_i2))%6;

    for(var _j = 0; _j < 2; _j++) {

      var _j2 = (2 + _j + (2*_i2))%6;

      var _solution = (-(
          _sliceNormal[_i2]*(_bbox[_i] - _sliceOrigin[_i2])
          +
          _sliceNormal[_i3]*(_bbox[_j2] - _sliceOrigin[_i3])
          )
          /
          _sliceNormal[_i4]
          )
          +
          _sliceOrigin[_i4]
          ;

      if((_solution >= _bbox[_j3] && _solution <= _bbox[_j3+1])
          ||
          (_solution <= _bbox[_j3] && _solution >= _bbox[_j3+1])) {

        var _sol = new Array();
        _sol[_i2] = _bbox[_i];
        _sol[_i3] = _bbox[_j2];
        _sol[_i4] = _solution;

        _solutionsIn.push(_sol);

      }
      else{

        var _sol = new Array();
        _sol[_i2] = _bbox[_i];
        _sol[_i3] = _bbox[_j2];
        _sol[_i4] = _solution;

        _solutionsOut.push(_sol);

      }
    }
  }

  return [_solutionsIn, _solutionsOut];
};

/**
 * Get XYToRAS transform and its inverse.
 *
 * @param {!Float32Array} _sliceNormal The slice normal.
 * @param {!Float32Array} _XYNormal The XY normal.
 *
 * @return The XY to RAS transform and its inverse.
 * @static
 */
parserTRK.xyrasTransform = function(_sliceNormal, _XYNormal){

  var _RASToXY = goog.vec.Mat4.createFloat32Identity();
    // no rotation needed if we are in the z plane already
  if(!goog.vec.Vec3.equals(_sliceNormal,_XYNormal)) {

    var _cp = _sliceNormal[2];
    var _teta = Math.acos(_cp);
    var _r = goog.vec.Vec3.createFloat32();
    goog.vec.Vec3.cross(_sliceNormal, _XYNormal, _r);
    goog.vec.Vec3.normalize(_r, _r);

    var a = Math.cos(_teta/2);
    var b = Math.sin(_teta/2)*_r[0];
    var c = Math.sin(_teta/2)*_r[1];
    var d = Math.sin(_teta/2)*_r[2];

    goog.vec.Mat4.setRowValues(_RASToXY,
        0,
        (a*a+b*b-c*c-d*d),
        2*(b*c-a*d),
        2*(b*d+a*c),
        0
        );
    goog.vec.Mat4.setRowValues(_RASToXY,
        1,
        2*(b*c+a*d),
        (a*a+c*c-b*b-d*d),
        2*(c*d-a*b),
        0
        );
    goog.vec.Mat4.setRowValues(_RASToXY,
        2,
        2*(b*d-a*c ),
        2*(c*d+a*b),
        (a*a+d*d-c*c-b*b),
        0
        );
    }


  var _XYToRAS = goog.vec.Mat4.createFloat32();
  goog.vec.Mat4.invert(_RASToXY, _XYToRAS);

  return [_RASToXY, _XYToRAS];
};

/**
 * Get bounding box given a point cloud.
 *
 * @param {!Array} _solutionsXY The slice origin in RAS space.
 *
 * @return The bounding box.
 * @static
 */
parserTRK.xyBBox = function(_solutionsXY){

  var _xyBBox = [Number.MAX_VALUE, -Number.MAX_VALUE,
   Number.MAX_VALUE, -Number.MAX_VALUE,
   Number.MAX_VALUE, -Number.MAX_VALUE];
  var i = 0;
  for (i = 0; i < _solutionsXY.length; ++i) {

    if(_solutionsXY[i][0] < _xyBBox[0]) {

      _xyBBox[0] = _solutionsXY[i][0];

    }

    if(_solutionsXY[i][0] > _xyBBox[1]) {

      _xyBBox[1] = _solutionsXY[i][0];

    }

    if(_solutionsXY[i][1] < _xyBBox[2]) {

      _xyBBox[2] = _solutionsXY[i][1];

    }

    if(_solutionsXY[i][1] > _xyBBox[3]) {

      _xyBBox[3] = _solutionsXY[i][1];

    }

    if(_solutionsXY[i][2] < _xyBBox[4]) {

      _xyBBox[4] = _solutionsXY[i][2];

    }

    if(_solutionsXY[i][2] > _xyBBox[5]) {

      _xyBBox[5] = _solutionsXY[i][2];

    }
  }

  return _xyBBox;
};
    
  // parse the header of the .TRK file
  // Documented here: http://trackvis.org/docs/?subsect=fileformat
  var header = {

    'id_string': parserTRK.scan('uchar', 6),
    'dim': parserTRK.scan('ushort', 3),
    'voxel_size': parserTRK.scan('float', 3),
    'origin': parserTRK.scan('float', 3),
    'n_scalars': parserTRK.scan('ushort'),
    'scalar_name': parserTRK.scan('uchar', 200),
    'n_properties': parserTRK.scan('ushort'),
    'property_name': parserTRK.scan('uchar', 200),
    'vox_to_ras': parserTRK.scan('float', 16),
    'reserved': parserTRK.scan('uchar', 444),
    'voxel_order': parserTRK.scan('uchar', 4),
    'pad2': parserTRK.scan('uchar', 4),
    'image_orientation_patient': parserTRK.scan('float', 6),
    'pad1': parserTRK.scan('uchar', 2),
    'invert_x': parserTRK.scan('uchar'),
    'invert_y': parserTRK.scan('uchar'),
    'invert_z': parserTRK.scan('uchar'),
    'swap_xy': parserTRK.scan('uchar'),
    'swap_yz': parserTRK.scan('uchar'),
    'swap_zx': parserTRK.scan('uchar'),
    'n_count': parserTRK.scan('uint'),
    'version': parserTRK.scan('uint'),
    'hdr_size': parserTRK.scan('uint')
  }; 

  //
  // parse the data

  // if n_count not provided, we parse the data until end of points
  /*var numberOfFibers = (header.n_count === 0) ? Infinity : header.n_count;
  var numberOfScalars = header.n_scalars;*/

  var m = new THREE.Matrix4()
  m.set(header.vox_to_ras[0],header.vox_to_ras[1],header.vox_to_ras[2],header.vox_to_ras[3],header.vox_to_ras[4],header.vox_to_ras[5],header.vox_to_ras[6],header.vox_to_ras[7],header.vox_to_ras[8],header.vox_to_ras[9],header.vox_to_ras[10],header.vox_to_ras[11],header.vox_to_ras[12],header.vox_to_ras[13],header.vox_to_ras[14],header.vox_to_ras[15]);
  
  // loop through all fibers
  var fibers = [];
  var lengths = [];
  var minLength = Infinity;
  var maxLength = -Infinity;

  var _numPoints = this.scan('uint', (this._data.byteLength - 1000) / 4);
  this.jumpTo(header.hdr_size);
  var _points = this.scan('float', (this._data.byteLength - 1000) / 4);

  var offset = 0;

  // keep track of the number of all points along all tracks
  var _totalPoints = 0;

  var i;
  var numberOfFibers = 0;
  for (i = 0; i < numberOfFibers; i++) {
    // if undefined, it means we have parsed all the data
    // (useful if n_count not defined or === 0)
    if(typeof(_numPoints[offset]) === 'undefined'){
      numberOfFibers = i;
      break;
    }

    var numPoints = _numPoints[offset];


    // console.log(numPoints, offset);


    var currentPoints = new THREE.Geometry()

    var length = 0.0;

    // loop through the points of this fiber
    for ( var j = 0; j < numPoints; j++) {

      // read coordinates
      var x = _points[offset + j * 3 + j * numberOfScalars + 1];
      var y = _points[offset + j * 3 + j * numberOfScalars + 2];
      var z =-1* _points[offset + j * 3 + j * numberOfScalars + 3];

      // console.log(x, y, z);

      // read scalars
      // var scalars = this.scan('float', header.n_scalars);

      // Convert coordinates to world space by dividing by spacing
      x = x / header.voxel_size[0];
      y = y / header.voxel_size[1];
      z = z / header.voxel_size[2];
    var vector=new THREE.Vector4( x,  y, z,1 )
    vector.applyMatrix4(m)
    vector.x-=0
    vector.y-=0
    vector.x*=1
    vector.y*=1
    currentPoints.vertices.push(new THREE.Vector3(vector.x, vector.y, vector.z) );

      // fiber length
      if (j > 0) {

        // if not the first point, calculate length
        var oldPoint = currentPoints.vertices[j - 1];
    var displacement=[Math.abs(vector.x - oldPoint.x), Math.abs(vector.y - oldPoint.y), Math.abs( vector.z- oldPoint.z)]
    curLength=Math.sqrt(displacement[0]*displacement[0] +
            displacement[1]*displacement[1] + displacement[2]*displacement[2]);
        length += curLength
    
    //adds in vertex color values
    if(j==1)
      currentPoints.colors.push( new THREE.Color( displacement[0]/curLength, displacement[1]/curLength, displacement[2]/curLength ))
    
    currentPoints.colors.push( new THREE.Color( displacement[0]/curLength, displacement[1]/curLength, displacement[2]/curLength ))
      }

      // increase the number of points if this is not the last track
      if (j < numPoints - 1) {
        _totalPoints += 6;
      }

    }
  currentPoints.computeBoundingBox();
  currentPoints.computeFaceNormals();
  currentPoints.computeVertexNormals();
    offset += numPoints * 3 + numPoints * numberOfScalars + 1;


    // read additional properties
    // var properties = this.scan('float', header.n_properties);

    // append this track to our fibers list
    fibers.push(currentPoints);
    // .. and also the length
    lengths.push(length);

  } // end of loop through all tracks

  

  // the scalar array
  var scalarArray = new Float32Array(_totalPoints);


  var _scalarIndex = 0;

  // now we have a list of fibers
  for (i = 0; i < numberOfFibers; i++) {

    // grab the current points of this fiber
    var points = fibers[i];
    var numberOfPoints = points.count;

    // grab the length of this fiber
    var length = lengths[i];

    minLength = Math.min(minLength, length);
    maxLength = Math.max(maxLength, length);

    for ( var j = 0; j < numberOfPoints - 1; j++) {

      // TODO min max check?

      // add the length (6 times since we added two points with each 3
      // coordinates)
      scalarArray[_scalarIndex++] = length;
      scalarArray[_scalarIndex++] = length;
      scalarArray[_scalarIndex++] = length;
      scalarArray[_scalarIndex++] = length;
      scalarArray[_scalarIndex++] = length;
      scalarArray[_scalarIndex++] = length;

    } // loop through points

  } // loop through fibers

  // make sure the vox_to_ras is not null. Else, set to identity
  var vox_to_ras_defined = false;
  for (i = 0; i < 16; i++) {
     if(header.vox_to_ras[i] != 0 ){
       vox_to_ras_defined = true;
       break;
     }
  }

  if(vox_to_ras_defined == false){
    header.vox_to_ras[0] = header.vox_to_ras[5] = header.vox_to_ras[10] = header.vox_to_ras[15] = 1;
  }



  // the object should be set up here, so let's fire a modified event
  object._points=fibers;
  console.log(fibers)
  modifiedEvent._object = object;
  modifiedEvent._container = container;
  modifiedEvent.marker=true;
  this.dispatchEvent(modifiedEvent);

};

var container = null;
var object = null;
var buffer = new ArrayBuffer(8);
var data = new Int32Array(buffer);
var flag = null;
var parser = parserTRK.call(parser, container, object, data, flag);

var prototype = parserTRK.prototype.parse.call(prototype, container, object, data, flag)
