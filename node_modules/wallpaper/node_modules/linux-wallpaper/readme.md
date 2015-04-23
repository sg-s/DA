# linux-wallpaper [![Build Status](https://travis-ci.org/sindresorhus/linux-wallpaper.svg?branch=master)](https://travis-ci.org/sindresorhus/linux-wallpaper)

> Get or set the desktop wallpaper on Linux


## Install

```
$ npm install --save linux-wallpaper
```


## Usage

```js
var linuxWallpaper = require('linux-wallpaper');

linuxWallpaper.set('unicorn.jpg', function (err) {
	console.log('done');
});

linuxWallpaper.get(function (err, imagePath) {
	console.log(imagePath);
	//=> '/home/sindresorhus/unicorn.jpg'
});
```


## API

### .get(callback)

#### callback(error, imagePath)

*Required*  
Type: `function`

##### imagePath

Type: `string`

Path to the current desktop wallpaper image.

### .set(imagePath, [callback])

##### imagePath

*Required*  
Type: `string`

Path to the image to set as the desktop wallpaper.


## CLI

See [`wallpaper`](https://github.com/sindresorhus/wallpaper#cli).


## License

MIT © [Sindre Sorhus](http://sindresorhus.com)
