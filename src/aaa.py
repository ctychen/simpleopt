import numpy as np
import plotly.graph_objects as go
import plotly.io as pio

# Create a 3D meshgrid using NumPy
x = np.linspace(-5, 5, 10)
y = np.linspace(-5, 5, 10)
z = np.linspace(-5, 5, 10)
X, Y, Z = np.meshgrid(x, y, z)

# Define the function to plot on the surface
F = np.sin(np.sqrt(X**2 + Y**2 + Z**2))

print(len(F))

print(F)

# print(X)

# print(Y)

# print(Z)

# Create a 3D surface plot with color mapped to function values
fig = go.Figure(data=go.Volume(
    x=X.flatten(),
    y=Y.flatten(),
    z=Z.flatten(),
    value=F.flatten(),
    isomin=np.min(F),
    isomax=np.max(F),
    opacity=0.2,
    surface_count=17,
    colorscale='Viridis'
))

# Set plot layout and axis labels
fig.update_layout(
    title='3D Surface Plot with Colormap',
    scene=dict(
        xaxis_title='X',
        yaxis_title='Y',
        zaxis_title='Z'
    )
)

# Show the plot
fig.show()



output_file = 'plottest.html'
pio.write_html(fig, output_file)




