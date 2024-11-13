def export_rea(filename, elements, boundaries):
    with open(filename, 'w') as file:
        n = elements.shape[0]
        # Add header
        with open('header.rea', 'r') as header_file:
            for line in header_file:
                file.write(line)

        # Mesh
        file.write('**MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n')
        file.write(f'{n:8d}   3   {n:8d} NEL NDIM NEL\n')
        for i in range(n):
            file.write(f'            ELEMENT{i+1:12d} [    1A]    GROUP     0\n')
            for j in range(3):
                file.write(f'{elements[i, 0, j]:15.7g}  {elements[i, 1, j]:15.7g}  '
                           f'{elements[i, 2, j]:15.7g}  {elements[i, 3, j]:15.7g}\n')
            for j in range(3):
                file.write(f'{elements[i, 4, j]:15.7g}  {elements[i, 5, j]:15.7g}  '
                           f'{elements[i, 6, j]:15.7g}  {elements[i, 7, j]:15.7g}\n')

        file.write('  ***** CURVED SIDE DATA *****\n')
        file.write('     0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n')
        file.write('  ***** BOUNDARY CONDITIONS ***** \n')
        file.write(' ***** FLUID   BOUNDARY CONDITIONS ***** \n')
        for i in range(n):
            for j in range(6):
                k = [item[2] for item in boundaries if item[0]==i and item[1] == j]
                if len(k) == 0:
                    k = 0
                else:
                    k = k[0]

                if k == 4:
                    file.write(f' W    {i+1:3d} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e}\n')
                elif k == 2:
                    file.write(f' v    {i+1:3d} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e}\n')
                elif k == 3:
                    file.write(f' O    {i+1:3d} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e}\n')
                else:
                    file.write(f' E    {i+1:3d} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e} {0:13.5e}\n')

        # Add footer
        with open('footer.rea', 'r') as footer_file:
            for line in footer_file:
                file.write(line)

