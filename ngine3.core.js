var Ngine3 = {}
//import Entity from 'ngine3.entity'
Ngine3.Components = {}
Ngine3.Systems = {}
Ngine3.Entity = function Entity () {
    this.id = (+new Date()).toString(16) + (Math.random() * 100000000 | 0).toString(16) + ECS.Entity.prototype._count
    Entity.prototype._count++
    this.components = {}
    return this
}

Ngine3.Entity.prototype._count = 0;

Ngine3.Entity.prototype.addComponent = function addComponent (component){
    this.components[component.name] = component
    return this
}

Ngine3.Entity.prototype.removeComponent = function removeComponent (componentName){
    var name = componentName
    if (typeof componentName === 'function') { 
        name = componentName.prototype.name
    }
    delete this.components[name]
    return this
}


Ngine3.Interface = function Interface(host){
    //alert(host.toString());
    var entity=this;
    var mouseXY;
    host.addEventListener('mousedown', mouseDown, false);
    host.addEventListener('mouseup', mouseUp, false);
    window.addEventListener('mousemove', mouseMove, false);
    host.addEventListener('click', click, false);
    window.addEventListener('keydown',keyDown,false);
    window.addEventListener('keyup',keyUp,false);
    window.addEventListener((/Firefox/i.test(navigator.userAgent))? "DOMMouseScroll" : "mousewheel",wheel,false);
    function click(e){
        entity.click(e);
    }
    function wheel(event){
        var delta = 0;
        if (!event) /* For IE. */
                event = window.event;
        if (event.wheelDelta) { /* IE/Opera. */
                delta = event.wheelDelta/120;
        } else if (event.detail) { /** Mozilla case. */
                /** In Mozilla, sign of delta is different than in IE.
                 * Also, delta is multiple of 3.
                 */
                delta = -event.detail/3;
        }
        /** If delta is nonzero, handle it.
         * Basically, delta is now positive if wheel was scrolled up,
         * and negative, if wheel was scrolled down.
         */
        if (delta)
                entity.wheel(event,delta);
    }
    function mouseDown(e){
        entity.mouseDown(e);
    }
    function mouseUp(e){
        entity.mouseUp(e);
    }
    function mouseMove(e){
        entity.mouseMove(e);
    }
    function keyDown(e){
        entity.keyDown(e);
    }
    function keyUp(e){
        entity.keyUp(e);
    }
    this.keyDown=function(e){}
    this.keyUp=function(e){}
    this.mouseDown=function(e){}
    this.mouseUp=function(e){}
    this.mouseMove=function(e){}
    this.click=function(e){}
    this.wheel=function(e,d){}
    this.getMouseXY=function(e){
        if(!e){
            return mouseXY;
        }
        var tempXY=getTempXY(e);
        if(dot(tempXY,tempXY)==0){
            return mouseXY;
        }
        else{
            var hostXY=findPos(host);
            var realXY=[
                tempXY[0]-hostXY[0],
                tempXY[1]-hostXY[1]
            ];
            mouseXY=realXY;
            return realXY;
        }
    };
	
	function getTempXY(e) {
	// Detect if the browser is IE or not.
	// If it is not IE, we assume that the browser is NS.
	  var IE = document.all?true:false
	  var tempX = 0
	  var tempY = 0

	  if (IE) { // grab the x-y pos.s if browser is IE
		tempX = event.clientX + document.body.scrollLeft
		tempY = event.clientY + document.body.scrollTop
	  } else {  // grab the x-y pos.s if browser is NS
		tempX = e.pageX
		tempY = e.pageY
	  }  
	  // catch possible negative values in NS4
	  if (tempX < 0){tempX = 0}
	  if (tempY < 0){tempY = 0}  
	  // show the position values in the form named Show
	  // in the text fields named MouseX and MouseY
	  //document.Show.MouseX.value = tempX
	  //document.Show.MouseY.value = tempY
	  return [tempX, tempY];
	}

	function findPos(obj) {
		var curleft = 0;
			var curtop = 0;
			if (obj.offsetParent) {
				do {
				curleft += obj.offsetLeft;
				curtop += obj.offsetTop;
				} while (obj = obj.offsetParent);
			}
			return [curleft,curtop];
	}
}