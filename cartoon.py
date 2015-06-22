import matplotlib.pyplot as plt
import matplotlib.path
import numpy as np

def get_axes_coordinates(text):
    ax = text.axes
    renderer = ax.figure.canvas.get_renderer()
    bbox = text.get_window_extent(renderer)
    axes_coordinates = bbox.inverse_transformed(ax.transAxes).get_points()
    return axes_coordinates

def get_axes_coordinates_from_list(text_list):
    (min_x, min_y), (max_x, max_y) = get_axes_coordinates(text_list[0])
    for text in text_list[1:]:
        (x0, y0), (x1, y1) = get_axes_coordinates(text)
        if x0 < min_x:
            min_x = x0
        if y0 < min_y:
            min_y = y0
        if x1 > max_x:
            max_x = x1
        if y1 > max_y:
            max_y = y1
    return (min_x, min_y), (max_x, max_y)

def bracket_below_text(text, color='black'):
    (x0, y0), (x1, y1) = get_axes_coordinates(text)

    bracket_x_buffer = 0.03
    bracket_height = 0.2
    bracket_y_buffer = -0.3

    text_height = y1 - y0
    text_width = x1 - x0
    left_edge = x0 - bracket_x_buffer * text_width
    right_edge = x1 + bracket_x_buffer * text_width
    bottom = y0 + bracket_y_buffer * text_height
    top = bottom + bracket_height * text_height

    path = [[left_edge, top],
            [left_edge, bottom],
            [right_edge, bottom],
            [right_edge, top],
           ]

    patch = plt.Polygon(path,
                        color=color,
                        fill=False,
                        closed=False,
                        clip_on=False,
                       )

    text.axes.add_patch(patch)

    x, y = text.get_position()
    
    info = {'x': x,
            'y': y,
            'center': x0 + 0.5 * text_width,
            'bottom':  bottom,
            'height': text_height,
            'width': text_width,
            'right_edge': right_edge,
            }
    
    return info

def rectangle_around_text(text, specific_nucleotide=None, **kwargs):
    (x0, y0), (x1, y1) = get_axes_coordinates(text)
    
    x_buffer = 0.03
    text_width = x1 - x0
    
    if specific_nucleotide != None:
        x_buffer = 0
        text_width /= 3
        
        x0 += text_width * specific_nucleotide
        x1 -= text_width * (2 - specific_nucleotide)
            
    text_height = y1 - y0
    rect_height = 1.1 * text_height
    rect_width = (1 + 2 * x_buffer) * text_width
    rect_x = x0 - x_buffer * text_width
     
    rect = plt.Rectangle((rect_x, y0), rect_width, rect_height, clip_on=False, **kwargs)
    
    text.axes.add_patch(rect)
    
    center = rect_x + 0.5 * rect_width
    top = y0 + rect_height
    
    return center, top, rect_height, rect_width

def rectangle_around_text_list(text_list, **kwargs):
    (x0, y0), (x1, y1) = get_axes_coordinates_from_list(text_list)
    
    x_buffer = 0.01
    text_width = x1 - x0
            
    text_height = y1 - y0
    rect_height = 1.1 * text_height
    rect_width = (1 + 2 * x_buffer) * text_width
    rect_x = x0 - x_buffer * text_width
     
    rect = plt.Rectangle((rect_x, y0), rect_width, rect_height, clip_on=False, **kwargs)
    
    text_list[0].axes.add_patch(rect)
    
    center = rect_x + 0.5 * rect_width
    top = y0 + rect_height
    
    return center, top, rect_height
    
def text_above_text(text, string, **kwargs):
    (x0, y0), (x1, y1) = get_axes_coordinates(text)
    above = text.axes.text((x0 + x1) / 2,
                           y1 + .2 * (y1 - y0),
                           string,
                           ha='center',
                           va='bottom',
                           **kwargs)
    return above

def connect_with_curve(ax, x0, x1, y, height, offset_text, font):
    points = [(x0, y),
              (0.8 * x0 + 0.2 * x1, y + 0.4 * height),   
              (0.2 * x0 + 0.8 * x1, y + 0.4 * height),
              (x1, y),
             ]

    codes = [matplotlib.path.Path.MOVETO,
             matplotlib.path.Path.CURVE4,
             matplotlib.path.Path.CURVE4,
             matplotlib.path.Path.CURVE4,
            ]

    path = matplotlib.path.Path(points, codes)
    patch = matplotlib.patches.PathPatch(path, facecolor='none', transform=ax.transAxes, clip_on='False')
    ax.add_patch(patch)
    
    font = font.copy()
    font['size'] = font['size'] * 0.75
    
    ax.text((x0 + x1) / 2,
            y + height,
            offset_text,
            ha='center',
            va='bottom',
            fontdict=font,
           )

def connect_with_arrow(ax, x0, x1, y, height, arrow_width, offset_text, font):
    points = [(x0, y),
              (x0, y + height),   
              (x1, y + height),
              (x1, y),
              (x1 - arrow_width, y + 0.5 * height),
              (x1, y),
              (x1 + arrow_width, y + 0.5 * height),
              (x1, y),
             ]

    codes = [matplotlib.path.Path.MOVETO,
             matplotlib.path.Path.LINETO,
             matplotlib.path.Path.LINETO,
             matplotlib.path.Path.LINETO,
             matplotlib.path.Path.MOVETO,
             matplotlib.path.Path.LINETO,
             matplotlib.path.Path.MOVETO,
             matplotlib.path.Path.LINETO,
            ]

    path = matplotlib.path.Path(points, codes)
    patch = matplotlib.patches.PathPatch(path, facecolor='none', transform=ax.transAxes, clip_on='False')
    ax.add_patch(patch)
    
    font = font.copy()
    font['size'] = font['size'] * 0.75
    font['family'] = 'serif'
    
    ax.text((x0 + x1) / 2,
            y + height,
            offset_text,
            ha='center',
            va='bottom',
            fontdict=font,
           )

class CodingSequence(object):
    def __init__(self, ax, buffered_codon_counts, gene_name,
                 font_size=40,
                 UTR_codons=0,
                 plot_up_to=30,
                 plot_end=False,
                ):
        self.ax = ax
        self.gene_name = gene_name
        self.font_size = font_size
        self.UTR_codons = UTR_codons
        self.plot_up_to = plot_up_to
        self.plot_end = plot_end
        
        cds_slice = slice(('start_codon', 0), ('stop_codon', 1))
        full_slice = slice(('start_codon', -self.UTR_codons), ('stop_codon', 1 + self.UTR_codons))

        self.codon_sequence = buffered_codon_counts[gene_name]['identities'][full_slice]
        self.codon_counts = buffered_codon_counts[gene_name]['relaxed'][full_slice]
        
        CDS_counts = buffered_codon_counts[gene_name]['relaxed'][full_slice]
        CDS_mean = np.mean(CDS_counts) 

        if plot_end:
            self.codon_seqeunce = self.codon_sequence[-plot_up_to:]
            self.codon_counts = self.codon_counts[-plot_up_to:]
        else:
            self.codon_sequence = self.codon_sequence[:plot_up_to]
            self.codon_counts = self.codon_counts[:plot_up_to]

        self.codon_enrichments = self.codon_counts / CDS_mean

    def draw_codon_sequence(self,
                            start_x, start_y,
                            highlight_offset=None,
                            specific_nucleotide=None,
                            highlight_sequence=None,
                            highlight_specific_position=None,
                            draw_brackets=True,
                            label=None,
                           ):
        
        sequence_positions = {-1: {'right_edge': None},
                             }
        texts = {}
        
        font = {'size': self.font_size,
                'family': 'monospace',
               }
        
        gap_between = 0.1
        
        if draw_brackets:
            bracket_color = 'black'
        else:
            bracket_color = 'white'
            
        for position, codon_id in enumerate(self.codon_sequence):
            previous_right_edge = sequence_positions[position - 1].get('right_edge')
            previous_width = sequence_positions[position - 1].get('width')
            
            if previous_right_edge == None:
                x = start_x
            else:
                x = previous_right_edge + gap_between * previous_width
            
            if self.plot_end:
                if len(self.codon_sequence) - position - 1 < self.UTR_codons:
                    alpha = 0.3
                else:
                    alpha = 1.0
            else:
                if position < self.UTR_codons:
                    alpha = 0.3
                else:
                    alpha = 1.0
            
            text = self.ax.text(x, start_y, codon_id, fontdict=font, alpha=alpha)        
            info = bracket_below_text(text, bracket_color)

            sequence_positions[position] = info
            texts[position] = text
            
        if highlight_offset != None:
            if specific_nucleotide != None:
                within_slice = slice(specific_nucleotide, specific_nucleotide + 1)
            else:
                within_slice = slice(None)

            positions_to_bold = set()

            color = 'black'
            if specific_nucleotide == None:
                if highlight_offset == 0:
                    color = 'red'
                elif highlight_offset == -1:
                    color = 'blue'
                elif highlight_offset == -2:
                    color = 'green'

            for position, codon_id in enumerate(self.codon_sequence):
                if not 0 < position + highlight_offset < len(self.codon_sequence):
                    continue
                    
                if highlight_specific_position != None and position != highlight_specific_position:
                    continue

                if self.codon_sequence[position + highlight_offset][within_slice] == highlight_sequence:
                    positions_to_bold.add(position)

                    offset_center, _, _, _ = rectangle_around_text(texts[position + highlight_offset],                                  
                                                                color=color,
                                                                alpha=0.4,
                                                                specific_nucleotide=specific_nucleotide,
                                                               )
                    position_center, top, height, width = rectangle_around_text(texts[position],
                                                                         color='black',
                                                                         alpha=1,
                                                                         fill=False,
                                                                        )
                    if highlight_offset != 0:
                        offset_text = '{0:+d}'.format(-highlight_offset)
                        connect_with_arrow(self.ax,
                                           offset_center,
                                           position_center,
                                           top,
                                           0.4 * height,
                                           0.15 * width,
                                           offset_text,
                                           font,
                                          )
                    
            sequence_positions['positions_to_bold'] = positions_to_bold
        else:
            sequence_positions['positions_to_bold'] = None
        
        first_codon = sequence_positions[0]
        name_x = first_codon['x'] - (1.5 if self.plot_end else 0.5) * first_codon['width']
        name_y = first_codon['y'] + 0.27 * first_codon['height']
         
        if label == None:
            label = '{0} CDS:'.format(self.gene_name)
        label_text = self.ax.text(name_x, name_y, label,
                                  size=self.font_size,
                                  ha='right',
                                  va='center',
                                  family='serif',
                                  color='blue',
                                 )
        coords = get_axes_coordinates(label_text)
        self.far_left = coords[0][0]
        
        if not self.plot_end or self.UTR_codons > 0 or hypothetical_sequence:
            last_codon = sequence_positions[len(self.codon_sequence) - 1]
            dots_x = last_codon['right_edge'] + (2 * gap_between) * last_codon['width']
            dots_y = last_codon['y'] 
            dots_text = self.ax.text(dots_x, dots_y, '...', fontdict=font)
            coords = get_axes_coordinates(dots_text)
            self.far_right = coords[1][0]
            sequence_positions['end dots'] = {'x': dots_x,
                                              'y': dots_y,
                                             }

        if self.plot_end:
            first_codon = sequence_positions[0]
            dots_x = first_codon['right_edge'] - (2 + gap_between) * first_codon['width']
            dots_y = first_codon['y'] 
            self.ax.text(dots_x, dots_y, '...', fontdict=font)
            sequence_positions['start dots'] = {'x': dots_x,
                                                'y': dots_y,
                                               }

        self.ax.axis('off')

        self.sequence_positions = sequence_positions
        # If counts aren't drawn, enrichments will use sequence_positions as count_positions
        self.count_positions = sequence_positions

    def draw_codon_counts(self, just_one_at=None):
        if just_one_at != None:
            codon_counts = np.zeros(self.plot_up_to, int)
            codon_counts[just_one_at] = 1
        else:
            codon_counts = self.codon_counts
        
        count_font_size = self.font_size * 0.5
        
        count_positions = {'positions_to_bold': self.sequence_positions['positions_to_bold'],
                          }
        
        gap_between = 0.1
        
        font = {'size': count_font_size}

        for i, count in enumerate(codon_counts):
            x = self.sequence_positions[i]['center']
            y = self.sequence_positions[i]['bottom'] - 0.9 * self.sequence_positions[i]['height']
            
            alpha = 1
            if count_positions['positions_to_bold'] != None and i not in count_positions['positions_to_bold']:
                alpha = 0.5
            
            text = self.ax.text(x, y, '{0:,}'.format(count), ha='center', va='top', fontdict=font, alpha=alpha)
            
            (x0, y0), (x1, y1) = get_axes_coordinates(text)
            
            count_positions[i] = {'center': self.sequence_positions[i]['center'],
                                  'bottom': y0,
                                  'height': y1 - y0,
                                 }
            
        name_x = self.sequence_positions[0]['x'] - 0.5 * self.sequence_positions[0]['width']
        name_y = count_positions[0]['bottom'] + 0.6 * count_positions[0]['height']
         
        self.ax.text(name_x, name_y, 'Read counts:', size=1.5 * count_font_size, family='serif', ha='right', va='center')
        
        if not self.plot_end or self.UTR_codons > 0:
            last_codon = self.sequence_positions[len(codon_counts) - 1]
            dots_x = self.sequence_positions['end dots']['x']
            dots_y = count_positions[0]['bottom'] + 0.1 * last_codon['height']
            self.ax.text(dots_x, dots_y, '...', size=count_font_size * 2, family='monospace')
            count_positions['end dots'] = {'x': dots_x,
                                           'y': dots_y,
                                          }
            
        self.count_positions = count_positions
            
    def draw_codon_enrichments(self):
        count_font_size = self.font_size * 0.5
        font = {'size': count_font_size}
        
        enrichment_positions = {}
        
        for i, enrichment in enumerate(self.codon_enrichments):
            x = self.count_positions[i]['center']
            y = self.count_positions[i]['bottom'] - 0.9 * self.count_positions[i]['height']
            
            enrichment_string = '{0:0.2f}'.format(enrichment)
            
            weight = 'normal'
            alpha = 1
            
            if self.count_positions['positions_to_bold'] != None:
                if i in self.count_positions['positions_to_bold']:
                    weight = 'bold'
                    alpha = 1
                else:
                    weight = 'normal'
                    alpha = 0.5
            
            text = self.ax.text(x, y, enrichment_string,
                                ha='center',
                                va='top',
                                weight=weight,
                                fontdict=font,
                                alpha=alpha,
                               )
            
            (x0, y0), (x1, y1) = get_axes_coordinates(text)
            
            enrichment_positions[i] = {'center': self.sequence_positions[i]['center'],
                                       'bottom': y0,
                                       'height': y1 - y0,
                                      }
            
        name_x = self.sequence_positions[0]['x'] - 0.5 * self.sequence_positions[0]['width']
        name_y = enrichment_positions[0]['bottom'] + 0.6 * enrichment_positions[0]['height']
         
        name_text = self.ax.text(name_x, name_y, 'Relative enrichment:',
                                 size=1.5 * count_font_size,
                                 family='serif',
                                 ha='right',
                                 va='center',
                                )
        coords = get_axes_coordinates(name_text)
        
        last_codon = self.sequence_positions[len(self.codon_counts) - 1]
        dots_x = self.count_positions['end dots']['x']
        dots_y = enrichment_positions[0]['bottom'] + 0.1 * last_codon['height']
        self.ax.text(dots_x, dots_y, '...', size=count_font_size * 2, family='monospace')

        self.bottom = coords[0][1]
            
    def draw_read(self,
                  position,
                  label_tRNA_sites='full',
                  label='Example read:',
                  y_gap=1.5,
                 ):
        font = {'size': self.font_size,
                'family': 'monospace',
               }
        
        texts = {}
        
        if self.plot_end:
            position = len(self.codon_sequence) - self.UTR_codons + position
        else:
            position = position + self.UTR_codons
        
        for i in range(position - 5, position + 5):
            codon_id = self.codon_sequence[i]
            if i == position + 4:
                codon_id = codon_id[:1]

            x = self.sequence_positions[i]['x']
            y = self.sequence_positions[i]['y'] + y_gap * self.sequence_positions[i]['height']
            
            texts[i] = self.ax.text(x, y, codon_id, color='white', fontdict=font)
            
        rectangle_around_text_list(texts.values(), color='black', fill=True, alpha=0.1)
        
        if label_tRNA_sites:
            tRNA_sites = [('A', 0, 'red'),
                          ('P', -1, 'blue'),
                          ('E', -2, 'green'),
                         ]
            
            for site, offset, color in tRNA_sites:
                rectangle_around_text(texts[position + offset], color='white', linewidth=0, fill=True)
                rectangle_around_text(texts[position + offset], color=color, linewidth=0, alpha=0.5, fill=True)
                if label_tRNA_sites == 'full':
                    above = text_above_text(texts[position + offset],
                                            '{0}-site'.format(site),
                                            color=color,
                                            size=self.font_size * 0.5,
                                           )
                    coords = get_axes_coordinates(above)
                    self.top = coords[1][1]
        
        for i in range(position - 5, position + 5):
            codon_id = self.codon_sequence[i]
            if i == position + 4:
                codon_id = codon_id[:1]

            x = self.sequence_positions[i]['x']
            y = self.sequence_positions[i]['y'] + y_gap * self.sequence_positions[i]['height']
            
            texts[i] = self.ax.text(x, y, codon_id, color='black', fontdict=font)
        
        rectangle_around_text_list(texts.values(), color='black', fill=False)

        (x0, y0), (x1, y1) = get_axes_coordinates(texts[position - 5])
        name_x = x0 - 0.2 * (x1 - x0)
        name_y = y0 + 0.5 * (y1 - y0)
         
        self.ax.text(name_x, name_y, label,
                     size=0.7 * self.font_size,
                     ha='right',
                     va='center',
                     family='serif',
                    )
        
    def rectangle_around_everything(self):
        rec = matplotlib.patches.Rectangle((self.far_left - 0.01, self.bottom - 0.04),
                                           self.far_right - self.far_left + 0.02,
                                           self.top - self.bottom + 0.08,
                                           fill=False,
                                           lw=1.2,
                                          )
        rec = self.ax.add_patch(rec)
        rec.set_clip_on(False)
