class Toolkit(object):
  
    def draw(self, plot, shapes):

        for name in shapes:
            shape = getattr(self, name, None)

            if shape is None:
                next

            shape(plot, shapes[name])

    def line(self, plot, lines, **kwargs):

        handle = plot.handle
        handle.line_color = kwargs.get('line_color','255 0 0')
        handle.line_width = kwargs.get('line_width','5')
        handle.line_style = kwargs.get('line_style','1')
        handle.clip       = str(kwargs.get('clip','1'))
        zorder            = kwargs.get('zorder','-1')

        for line in lines:

            if isinstance(lines, dict):

                collection = dict(kwargs)
                collection.update(lines[line])
                self.flight_path(plot, collection['data'], **collection)

            else:

                handle.line = line

                plot.cmd("""
                  set rgb $* $line_color
                  set line $* $line_style $line_width
                  set CLIP $clip
                  draw line $line
                """, zorder=zorder
                )

    def polygon(self, plot, polygons, **kwargs):

        handle = plot.handle
        handle.line_color = kwargs.get('line_color','255 0 0')
        handle.line_width = kwargs.get('line_width','5')
        handle.line_style = kwargs.get('line_style','1')
        zorder            = kwargs.get('zorder','-1')

        for poly in polygons:

            if isinstance(polygons, dict):

                collection = dict(kwargs)
                collection.update(polygons[poly])
                self.polygon(plot, collection['data'], **collection)

            else:

                handle.poly = poly

                plot.cmd("""
                  set rgb $* $line_color
                  set line $* $line_style $line_width
                  draw polyf $poly
                """, zorder=zorder
                )

    def rectangle(self, plot, rectangles, **kwargs):

        handle = plot.handle
        handle.line_color = kwargs.get('line_color','255 0 0')
        handle.fill_color = kwargs.get('fill_color',None)
        handle.line_width = kwargs.get('line_width','5')
        handle.line_style = kwargs.get('line_style','1')
        zorder            = kwargs.get('zorder','-1')

        for rect in rectangles:

            if isinstance(rect, dict):

                collection = dict(kwargs)
                collection.update(rect)
                self.rectangle(plot, collection['data'], **collection)

            elif isinstance(rectangles, dict):

                collection = dict(kwargs)
                collection.update(rectangles[rect])
                self.rectangle(plot, collection['data'], **collection)

            else:

                handle.rect = rect

                if handle.fill_color:
                    plot.cmd("""
                      set rgb $* $fill_color
                      set line $* $line_style $line_width
                      draw recf $rect
                    """, zorder=zorder
                    )

                if handle.line_color:
                    plot.cmd("""
                      set rgb $* $line_color
                      set line $* $line_style $line_width
                      draw rec $rect
                    """, zorder=zorder
                    )

    def string(self, plot, strings, **kwargs):

        handle = plot.handle
        handle.line_color = kwargs.get('line_color','255 255 255')
        handle.line_width = kwargs.get('line_width','5')
        handle.str_size   = kwargs.get('str_size','0.1 0.1')
        handle.position   = kwargs.get('position','c')
        handle.rotation   = kwargs.get('rotation','0')
        handle.font       = '$' + kwargs.get('font','regular')
        handle.clip       = str(kwargs.get('clip','1'))
        zorder            = kwargs.get('zorder','-1')

        for s in strings:

            if isinstance(s, dict):

                collection = dict(kwargs)
                collection.update(s)
                self.string(plot, collection['data'], **collection)

            elif isinstance(strings, dict):

                collection = dict(kwargs)
                collection.update(strings[s])
                self.string(plot, collection['data'], **collection)

            else:

                handle.str = s

                plot.cmd("""
                  set rgb $* $line_color
                  set font $font
                  set CLIP $clip
                  set string $* $position $line_width $rotation
                  set strsiz $str_size
                  draw string $str
                  """, zorder=zorder
                  )

    def mark(self, plot, marks, **kwargs):

        handle = plot.handle
        handle.fill_color   = kwargs.get('fill_color','255 0 0')
        handle.line_color   = kwargs.get('line_color',handle.fill_color)
        handle.line_width   = kwargs.get('line_width','5')
        handle.line_style   = kwargs.get('line_style','1')
        handle.size         = kwargs.get('size','.15')
        handle.type         = kwargs.get('mark','3')
        zorder              = kwargs.get('zorder','-1')

        for mark in marks:

            if isinstance(marks, dict):

                collection = dict(kwargs)
                collection.update(marks[mark])
                self.mark(plot, collection['data'], **collection)

            else:

                handle.mark = mark

                plot.cmd("""
                  set rgb $* $line_color
                  set line $* $line_style $line_width
                  draw mark $type $mark $size
                """, zorder=zorder
                )

    def station_mark(self, plot, marks, **kwargs):

        handle = plot.handle
        handle.line_color   = kwargs.get('line_color','255 0 0')
        handle.fill_color   = kwargs.get('fill_color','255 255 255')
        handle.line_width   = kwargs.get('line_width','5')
        handle.line_style   = kwargs.get('line_style','1')
        handle.inner_size   = kwargs.get('inner_size','.15')
        handle.outer_size   = kwargs.get('outer_size','.20')
        handle.outer_line   = kwargs.get('outer_line','.25')
        zorder              = kwargs.get('zorder','-1')

        if not handle.outer_line:
            handle.outer_line = '--auto'

        mark_type           = kwargs.get('mark_type','3 2')
        mark_type           = [c for c in mark_type if c != ' ']

        handle.mark_type1 = mark_type[0]
        handle.mark_type2 = mark_type[1]

        for mark in marks:

            if isinstance(marks, dict):

                collection = dict(kwargs)
                collection.update(marks[mark])
                self.station_mark(plot, collection['data'], **collection)

            else:

                mark = [s for s in mark.split() if s != ' ']

                handle.navigate = ''
                handle.mark     = ' '.join(mark[0:2])
                if len(mark) > 2: handle.navigate = ' '.join(mark[2:])

                plot.cmd("""
                  set rgb $* $fill_color
                  set line $* $line_style $line_width
                  draw mark $mark_type1 $mark $outer_size
                  set rgb $* $line_color
                  set line $* $line_style $line_width
                  draw mark $mark_type2 $mark $outer_line
                  draw mark $mark_type2 $mark $inner_size $navigate
                """, zorder=zorder
                )

    __call__ = draw

    flight_path = line
