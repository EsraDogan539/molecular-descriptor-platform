# Embedding ChalMolDB in another website

To show ChalMolDB inside another page (for example a department website), use an iframe with both parameters
below. `embed=true` hides Streamlit's own toolbar and padding; `embedded=1` lets ChalMolDB keep that mode when
visitors move between its pages.

```html
<iframe src="https://chalmoldb.streamlit.app/?embed=true&embedded=1"
        style="width:100%; height:900px; border:0;" title="ChalMolDB"></iframe>
<p><a href="https://chalmoldb.streamlit.app" target="_blank" rel="noopener">Open ChalMolDB in a new tab</a></p>
```

Notes:
- Keep the "open in a new tab" link: the hosted app sleeps after a period without visitors and needs up to a
  minute to wake, and some browsers restrict embedded pages.
- The credit line (developers, affiliation, how to cite) is part of the app and is shown in the embedded view.
- Streamlit documentation: https://docs.streamlit.io/deploy/streamlit-community-cloud/share-your-app/embed-your-app
